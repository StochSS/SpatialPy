# pylint: disable-msg=C0301
# pylint: disable-msg=C0103

import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
import warnings
import uuid


import numpy
import scipy.io
import scipy.sparse

from model import *

import inspect

try:
    # This is only needed if we are running in an Ipython Notebook
    import IPython.display
except:
    pass

try:
    import h5py
except:
    raise Exception("PyURDME requires h5py.")

try:
    import dolfin
    import mshr
except:
    raise Exception("PyURDME requires FeniCS/Dolfin.")

try:
    dolfin.parameters["linear_algebra_backend"] = "uBLAS"
except:
    dolfin.parameters["linear_algebra_backend"] = "Eigen"

import pickle
import json
import functools

# module-level variable to for javascript export in IPython/Jupyter notebooks
__pyurdme_javascript_libraries_loaded = False
def load_pyurdme_javascript_libraries():
    global __pyurdme_javascript_libraries_loaded
    if not __pyurdme_javascript_libraries_loaded:
        __pyurdme_javascript_libraries_loaded = True
        import os.path
        import IPython.display
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/three.js')) as fd:
            bufa = fd.read()
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/render.js')) as fd:
            bufb = fd.read()
        with open(os.path.join(os.path.dirname(__file__),'data/three.js_templates/js/OrbitControls.js')) as fd:
            bufc = fd.read()
        IPython.display.display(IPython.display.HTML('<script>'+bufa+bufc+bufb+'</script>'))


def deprecated(func):
    '''This is a decorator which can be used to mark functions
     as deprecated. It will result in a warning being emitted
     when the function is used.'''

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
             "Call to deprecated function {}.".format(func.__name__),
             category=DeprecationWarning,
             filename=func.func_code.co_filename,
             lineno=func.func_code.co_firstlineno + 1
         )
        return func(*args, **kwargs)
    return new_func


# Set log level to report only errors or worse
dolfin.set_log_level(dolfin.ERROR)
import logging
logging.getLogger('FFC').setLevel(logging.ERROR)
logging.getLogger('UFL').setLevel(logging.ERROR)

class URDMEModel(Model):
    """
        An URDME Model extends Model with spatial information and methods to create URDME solver input.
        TODO: Documentiation.
    """

    def __init__(self, name=""):
        Model.__init__(self, name)

        # Currently not used
        self.geometry = None
        #
        self.sd = []
        self.sd_initialized = False

        self.mesh = None
        self.xmesh = None
        self.stiffness_matrices = None
        self.mass_matrices = None

        # subdomains is a list of MeshFunctions with subdomain marker information
        self.subdomains = OrderedDict()
        self.old_style_subdomain = False

        # This dictionary hold information about the subdomains each species is active on
        self.species_to_subdomains = {}
        self.tspan = None

        # URDMEDataFunction objects to construct the data vector.
        self.listOfDataFunctions = []

        # Volume of each voxel in the dolfin dof ordering (not vertex ordering).
        self.dofvol = None

    def __getstate__(self):
        """ Used by pickle to get state when pickling. Because we
            have Swig wrappers to extension modules, we need to remove some instance variables
            for the object to pickle. """

        # Filter out any instance variable that is not picklable...
        state = {}
        for key, item in self.__dict__.items():
            if key == "subdomains":
                sddict = OrderedDict()
                for sdkey, sd_func in item.items():
                    tmpfile = tempfile.NamedTemporaryFile(suffix=".xml")
                    dolfin.File(tmpfile.name) << sd_func
                    tmpfile.seek(0)
                    sddict[sdkey] = tmpfile.read()
                    tmpfile.close()

                state[key] = sddict
            elif key in ["stiffness_matrices", "mass_matrices", "xmesh"]:
                state[key] = None
            else:
                state[key] = item


        return state

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """

        self.__dict__ = state

        if 'subdomains' in state:
            # Recreate the subdomain functions
            try:
                sddict = OrderedDict()
                for sdkey, sd_func_str in state["subdomains"].items():
                    fd = tempfile.NamedTemporaryFile(suffix=".xml")
                    fdname = fd.name
                    fd.write(sd_func_str)
                    fd.seek(0)
                    fd_in = dolfin.File(fdname)
                    func = dolfin.MeshFunction("size_t", self.__dict__["mesh"])
                    fd_in >> func
                    sddict[sdkey] = func
                    fd.close()
                self.__dict__["subdomains"] = sddict
            except Exception as e:
                raise Exception("Error unpickling model, could not recreate the subdomain functions"+str(e))

        self.create_extended_mesh()

    def add_data_function(self, data_function):
        """ Add a URDMEDataFunction object to this object. """
        if isinstance(data_function, URDMEDataFunction):
            self.listOfDataFunctions.append(data_function)
        else:
            raise Exception("data_function not of type URDMEDataFunction")

    def __initialize_species_map(self):
        i = 0
        self.species_map = {}
        for S in self.listOfSpecies:
            self.species_map[S] = i
            i = i + 1

    def get_species_map(self):
        """ Get the species map, name to index. """
        if not hasattr(self, 'species_map'):
            self.__initialize_species_map()

        return self.species_map

    def add_subdomain(self, subdomain, subdomain_id=None):
        """ Add a subdomain definition to the model.  By default, all regions are set to
        subdomain 1.
        Args:
            subdomain: an instance of a 'dolfin.SubDomain' subclass.
            id: an int, the identifier for this subdomain.
        """
        if subdomain_id is None and isinstance(subdomain, dolfin.cpp.mesh.MeshFunctionSizet):
            # Old style, for backwards compatability
            if not subdomain.dim() in self.subdomains.keys():
                self.subdomains[subdomain.dim()] = subdomain
                self.old_style_subdomain = True
            else:
                raise ModelException("Failed to add subdomain function of dim "+str(subdomain.dim())+". Only one subdomain function of a given dimension is allowed.")
        else:
            # New style
            if self.old_style_subdomain:
                raise ModelException("Can't mix old and new style subdomains")
            if not issubclass(subdomain.__class__, dolfin.SubDomain):
                raise ModelException("'subdomain' argument to add_subdomain() must be a subclass of dolfin.SubDomain")
            if subdomain_id is None or not isinstance(subdomain_id, int):
                raise ModelException("'id' argument to add_subdomain() must be an int")
            self.subdomains[subdomain_id] = subdomain


    ################################

    def _subdomains_to_threejs(self, subdomains={1:'blue'}):
        """ Export threejs code to plot the mesh with edges colored for the listed subdomains.
            Input is a dictionary with subdomain index:color pairs, output is a single json three.js mesh
            with the subdomains colored according to the input colors. """
        sd = self.get_subdomain_vector()
        c = ['black']*len(sd)

        for i, s in enumerate(sd):
            try:
                c[i] = subdomains[int(s)]
            except KeyError:
                pass

        jsondoc = self.mesh.export_to_three_js(colors = c)
        return jsondoc

    def _subdomains_to_html(self, filename, sd=1):
        sd = self.get_subdomain_vector()
        c = _compute_colors(sd)
        self.mesh._ipython_display_(filename, colors=c)
    
    def write_stochss_subdomain_file(self, filename="stochss_sdfile.txt"):
        # Write mesh and subdomain files for the StochSS UI
        sd = self.get_subdomain_vector()
        with open(filename,'w') as fd:
            for ndx, val in enumerate(sd):
                fd.write("{0},{1}\n".format(ndx, val))

    def display_mesh(self, subdomains, width=500, height=375, camera=[0,0,1]):
        ''' WebGL display of the wireframe mesh.'''
        load_pyurdme_javascript_libraries()
        if isinstance(subdomains, int):
            jstr = self._subdomains_to_threejs(subdomains={1:'blue', subdomains:'red'})
        elif isinstance(subdomains, list):
            sd_in = {1:'blue'}
            for i in subdomains:
                sd_in[i] = 'red'
            jstr = self._subdomains_to_threejs(subdomains=sd_in)
        elif isinstance(subdomains, dict):
            jstr = self._subdomains_to_threejs(subdomains=subdomains)
        hstr = None
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/mesh.html", 'r') as fd:
            hstr = fd.read()
        if hstr is None:
            raise Exception("could note open template mesh.html")
        hstr = hstr.replace('###PYURDME_MESH_JSON###', jstr)
        # Create a random id for the display div. This is to avioid multiple plots ending up in the same
        # div in Ipython notebook
        displayareaid = str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###', displayareaid)
        # ###CAMERA_X###, ###CAMERA_Y###, ###CAMERA_Z###
        hstr = hstr.replace('###CAMERA_X###',str(camera[0]))
        hstr = hstr.replace('###CAMERA_Y###',str(camera[1]))
        hstr = hstr.replace('###CAMERA_Z###',str(camera[2]))
        html = '<div style="width: {0}px; height: {1}px;" id="{2}" ></div>'.format(width, height, displayareaid)
        IPython.display.display(IPython.display.HTML(html+hstr))



    ################################

    def create_stoichiometric_matrix(self):
        """ Generate a stoichiometric matrix in sparse CSC format. """

        if not hasattr(self, 'species_map'):
            self.__initialize_species_map()
        if self.get_num_reactions() > 0:
            ND = numpy.zeros((self.get_num_species(), self.get_num_reactions()))
            for i, r in enumerate(self.listOfReactions):
                R = self.listOfReactions[r]
                reactants = R.reactants
                products  = R.products

                for s in reactants:
                    ND[self.species_map[s], i] -= reactants[s]
                for s in products:
                    ND[self.species_map[s], i] += products[s]

            N = scipy.sparse.csc_matrix(ND)
        else:
            N = numpy.zeros((self.get_num_species(), self.get_num_reactions()))

        return N

    def create_dependency_graph(self):
        """ Construct the sparse dependency graph. """
        # We cannot safely generate a dependency graph (without attempting to analyze the propensity string itself)
        # if the model contains custom propensities.
        mass_action_model = True
        for name, reaction in self.listOfReactions.items():
            if not reaction.massaction:
                GF = numpy.ones((self.get_num_reactions(), self.get_num_reactions() + self.get_num_species()))
                mass_action_model = False

        if mass_action_model:
            GF = numpy.zeros((self.get_num_reactions(), self.get_num_reactions() + self.get_num_species()))
            species_map = self.get_species_map()

            involved_species = []
            reactants = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                temp2 = []
                for s in reaction.reactants:
                    temp.append(species_map[s])
                    temp2.append(species_map[s])
                for s in reaction.products:
                    temp.append(species_map[s])
                involved_species.append(temp)
                reactants.append(temp2)

            species_to_reactions = []
            for species in self.listOfSpecies:
                temp = []
                for j, x in enumerate(reactants):
                    if species_map[species] in x:
                        temp.append(j)
                species_to_reactions.append(temp)


            reaction_to_reaction = []
            for name, reaction in self.listOfReactions.items():
                temp = []
                for s in reaction.reactants:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+species_to_reactions[species_map[s]]

                for s in reaction.products:
                    if species_to_reactions[species_map[s]] not in temp:
                        temp = temp+ species_to_reactions[species_map[s]]

                temp = list(set(temp))
                reaction_to_reaction.append(temp)

            # Populate G
            for j, spec in enumerate(species_to_reactions):
                for s in spec:
                    GF[s, j] = 1

            for i,reac in enumerate(reaction_to_reaction):
                for r in reac:
                    GF[r, self.get_num_species()+i] = 1


        try:
            G = scipy.sparse.csc_matrix(GF)
        except Exception as e:
            G = GF

        return G


    def timespan(self, tspan):
        """ Set the time span of simulation. """
        self.tspan = tspan


    def _initialize_default_subdomain(self):
        """" Create a default subdomain function. The default is all voxels belong
             to subdomain 1.
        """

        subdomain = dolfin.MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1)
        subdomain.set_all(1)
        self.add_subdomain(subdomain)

    def _initialize_species_to_subdomains(self):
        """ Initialize the species mapping to subdomains. The default
            is that a species is active in all the defined subdomains.
        """

        sds = list(numpy.unique(self.get_subdomain_vector()))
        # This conversion is necessary for UFL not to choke on the subdomain ids.
        for i, sd in enumerate(sds):
            sds[i] = int(sd)
        try:
            sds.remove(0)
        except ValueError:
            pass

        # If a species is not present as key in the species_to_subdomain mapping,
        # we label it as active in all subdomains
        for spec_name in self.listOfSpecies:
            species = self.listOfSpecies[spec_name]
            if species not in self.species_to_subdomains.keys():
                self.species_to_subdomains[species] = sds


    def restrict(self, species, subdomains):
        """ Restrict the diffusion of a species to a subdomain. """
        self.species_to_subdomains[species] = subdomains


    def set_subdomain_vector(self, sd):
        """ Explicitly set the subdomain vector from an array. """
        self.sd = sd
        self.sd_initialized = True

    def get_subdomain_vector(self, subdomains=None):
        """ Create the 'sd' vector. 'subdomains' is a dolfin FacetFunction,
            and if no subdomain input is specified, they voxels default to
            subdomain 1. """
        if self.sd_initialized:
            return self.sd

        # We need to make sure that the highest dimension is applied
        # first, otherwise the cell level will overwrite all markings
        # applied on boundaries.
        if not hasattr(self,'xmesh'):
            self.create_extended_mesh()

        self.mesh.init()

        sd = numpy.ones(self.mesh.get_num_voxels())

        if len(self.subdomains) == 0:
            self.sd = sd
        else:
            if self.old_style_subdomain:
                subdomains = self.subdomains
            else:
                subdomains = OrderedDict()
                sdvec = dolfin.MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1)
                sdvec.set_all(1)
                for id, inst in self.subdomains.iteritems():
                    inst.mark(sdvec, id)
                subdomains[sdvec.dim()] = sdvec

            for dim, subdomain in subdomains.items():
                if dim == 0:
                    # If we define subdomains on vertices, ONLY use those.
                    # Then it is a direct copy to the sd
                    for ndx, val in enumerate(subdomain):
                        sd[ndx] = val
                    break
                else:
                    # Map all facet labels to vertex labels
                    tovertex = self.mesh.topology()(dim, 0)
                    for i in range(subdomain.size()):
                        for vtx in tovertex(i):
                            if subdomain[i] != 0: # TODO: Temporary hack to fix issue with Gmesh facet_region files.
                                sd[vtx] = subdomain[i]

        self.sd = sd
        self.sd_initialized = True
        return self.sd

    def initialize_initial_condition(self):
        """ Create all-zeros inital condition matrix. """
        
        ns = self.get_num_species()
        if self.xmesh == None:
            self.create_extended_mesh()
        nv = self.mesh.get_num_voxels()
        self.u0 = numpy.zeros((ns, nv))

    def create_extended_mesh(self):
        """ Extend the primary mesh with information about the degrees of freedom. """

        xmesh = URDMEXmesh()
        # Construct a species map (dict mapping model species name to an integer index)
        species_map = self.get_species_map()
        # Initialize the function spaces and dof maps.
        for spec in self.listOfSpecies:
            species = self.listOfSpecies[spec]
            spec_name = species.name
            spec_index = species_map[spec_name]
            xmesh.function_space[spec_name] = self.mesh.get_function_space()
            xmesh.vertex_to_dof_map[spec_name] = dolfin.vertex_to_dof_map(xmesh.function_space[spec_name])
            xmesh.vertex_to_dof_map[spec_name] = len(self.listOfSpecies) * xmesh.vertex_to_dof_map[spec_name] + spec_index
            xmesh.dof_to_vertex_map[spec_name] = dolfin.dof_to_vertex_map(xmesh.function_space[spec_name])


        xmesh.vertex = self.mesh.coordinates()
        self.xmesh = xmesh


    # Some utility routines to set initial conditions
    def set_initial_condition_scatter(self, spec_init, subdomains=None):
        """ Scatter an initial number of molecules over the voxels in a subdomain. """

        if not hasattr(self,"u0"):
            self.initialize_initial_condition()

        if not hasattr(self, 'xmesh'):
            self.create_extended_mesh()

        self._initialize_species_to_subdomains()

        self.get_subdomain_vector()

        for species in spec_init:

            if subdomains is None:
                subdomains = self.species_to_subdomains[species]

            spec_name = species.name
            num_spec = spec_init[species]
            species_map = self.get_species_map()
            specindx = species_map[spec_name]

            sd = self.sd
            table = []
            for i, ind in enumerate(sd):
                if ind in subdomains:
                    table.append(i)

            ltab = len(table)
            if ltab == 0:
                raise ModelException("set_initial_condition_scatter: No voxel in the given subdomains "+str(subdomains)+", check subdomain marking.")

            for mol in range(num_spec):
                vtx = numpy.random.randint(0, ltab)
                ind = table[vtx]
                self.u0[specindx, ind] += 1

    def set_initial_condition_distribute_uniformly(self, spec_init, subdomains=None):
        """ Place the same number of molecules of the species in each voxel. """
        if not hasattr(self, "u0"):
            self.initialize_initial_condition()

        if not hasattr(self, 'xmesh'):
            self.create_extended_mesh()

        self._initialize_species_to_subdomains()

        species_map = self.get_species_map()
        for spec in spec_init:
            if subdomains is None:
                subdomains = self.species_to_subdomains[spec]
            spec_name = spec.name
            num_spec = spec_init[spec]
            specindx = species_map[spec_name]
            for ndx in range(len(self.sd)):
                if self.sd[ndx] in subdomains:
                    self.u0[specindx, ndx] = num_spec
    
    
    def set_initial_condition_place_near(self, spec_init, point=None, add=False):
        """ Place all molecules of kind species in the voxel nearest a given point. The species existing previously in this voxel are reset if add is set to False"""


        if not hasattr(self, "u0"):
            self.initialize_initial_condition()

        if not hasattr(self, 'xmesh'):
            self.create_extended_mesh()

        for spec in spec_init:
            spec_name = spec.name
            num_spec = spec_init[spec]

            # Find the voxel with center (vertex) nearest to the point
            ix = self.mesh.closest_vertex(point)
            species_map = self.get_species_map()
            specindx = species_map[spec_name]
            self.u0[specindx, ix] = (self.u0[specindx, ix] if add else 0) + num_spec

    def set_initial_condition_place_voxel(self, spec_init, voxel,add=False):
        """Place all molecules of kind species in the given voxel. The species existing previously in this voxel are reset if add is set to False"""

        if not hasattr(self, "u0"):
            self.initialize_initial_condition()

        if not hasattr(self, 'xmesh'):
            self.create_extended_mesh()

        for spec in spec_init:
            spec_name = spec.name
            num_spec = spec_init[spec]

            species_map = self.get_species_map()
            specindx = species_map[spec_name]
            self.u0[specindx, voxel] = (self.u0[specindx,voxel] if add else 0) + num_spec

    def create_system_matrix(self):
        """ Create the system (diffusion) matrix for input to the URDME solvers. The matrix
            is built by concatenating the individually assembled matrices for each of the species,
            and multiplying with the lumped mass matrix (which define the volume of the voxels).
            Negative off-diagonal elements in the matrix are set to zero, and the diagonal is renormalized
            in order to assure that the returned matrix is a Markov transition matrix.
            Returns a dictionary containing the volumes of the subvolumes, the system diffusion matrix
            and the fraction of the mass of the negative off-diagonal elements that has been filtered out.
            """

        import time

        # Check if the individual stiffness and mass matrices (per species) have been assembled, otherwise assemble them.
        if self.stiffness_matrices is not None and self.mass_matrices is not None:
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices
        else:
            if self.mesh is None:
                raise ModelException("This model has no mesh, can not create system matrix.")
            matrices = self.assemble()
            self.stiffness_matrices = matrices['K']
            self.mass_matrices = matrices['M']
            stiffness_matrices = self.stiffness_matrices
            mass_matrices = self.mass_matrices

        # Make a dok matrix of dimension (Ndofs,Ndofs) for easier manipulation
        i = 1
        Mspecies = len(self.listOfSpecies)
        if Mspecies == 0:
            raise ModelException("The model has no species, can not create system matrix.")
        
        # Use dolfin 'dof' number of voxels, not the number of verticies
        Nvoxels = self.mesh.get_num_dof_voxels()
        Ndofs = Nvoxels*Mspecies

        # Create the volume vector by lumping the mass matrices
        vol = numpy.zeros((Ndofs, 1))
        spec = 0

        xmesh = self.xmesh

        for species, M in mass_matrices.iteritems():

            rows, cols, vals = dolfin.as_backend_type(M).data()
            SM = scipy.sparse.csr_matrix((vals, cols, rows))
            vols = SM.sum(axis=1)

            spec = self.species_map[species]
            for j in range(len(vols)):
                vx = j
                dof = Mspecies*vx+spec
                vol[dof, 0] = vols[j]

        # This is necessary in order for the array to have the right dimension (Ndofs,1)
        vol = vol.flatten()

        # Assemble one big matrix from the individual stiffness matrices. Multiply by the inverse of
        # the lumped mass matrix, filter out any entries with the wrong sign and renormalize the columns.
        spec = 0
        positive_mass = 0.0
        total_mass = 0.0


        sd = self.get_subdomain_vector()
        sd_vec_dof = numpy.zeros(self.mesh.get_num_dof_voxels())
        vertex_to_dof = dolfin.vertex_to_dof_map(self.mesh.get_function_space())
        for ndx, sd_val in enumerate(sd):
            sd_vec_dof[vertex_to_dof[ndx]] = sd_val
        sd = sd_vec_dof

        tic  = time.time()
        # If a volume is zero, we need to set it to 1.
        vi = vol+(vol<=0.0)

        S = scipy.sparse.dok_matrix((Ndofs, Ndofs))


        keys = []
        values = []

        for species, K in stiffness_matrices.iteritems():

            rows, cols, vals = dolfin.as_backend_type(K).data()

            # Filter the matrix: get rid of all elements < 0 (inlcuding the diagonal)
            vals *= vals < 0
            Kcrs = scipy.sparse.csr_matrix((vals, cols, rows))

            sdmap  = self.species_to_subdomains[self.listOfSpecies[species]]

            # Filter the matrix: get rid of all elements < 0 (inlcuding the diagonal)
            Kdok = Kcrs.todok()


            for ind, val in Kdok.iteritems():

                ir = ind[0]
                ij = ind[1]

                # Check if this is an edge that the species should diffuse along,
                # if not, set the diffusion coefficient along this edge to zero. This is
                # equivalent to how boundary species are handled in legacy URDME.
                if sd[ir] not in sdmap:
                    val = 0.0

                S[Mspecies*ir+spec, Mspecies*ij+spec] = -val/vi[Mspecies*ij+spec]

            spec = spec + 1

        sumcol = S.tocsr().sum(axis=0)
        S.setdiag(-numpy.array(sumcol).flatten())

        D = S.tocsc()

        if total_mass == 0.0:
            return {'vol':vol, 'D':D, 'relative_positive_mass':None}
        else:
            return {'vol':vol, 'D':D, 'relative_positive_mass':positive_mass/total_mass}


    def validate(self, urdme_solver_data):
        """ Validate the model data structures.
            validate should be called prior to writing the model to the solver input file.
            The core solvers do very limited error checking of the input.
        """

        for spec_name, species in self.listOfSpecies.items():
            if 0 in self.species_to_subdomains[species]:
                raise ModelException("Subdomain number 0 is reserved. Please check your model.")

        # Check that all the columns of the system matrix sums to zero (or close to zero). If not, it does
        # not define a Markov process and the solvers might segfault or produce erraneous results.
        colsum = numpy.abs(urdme_solver_data['D'].sum(axis=0))
        colsum = colsum.flatten()
        maxcolsum = numpy.argmax(colsum)
        if colsum[0, maxcolsum] > 1e-10:
            D = urdme_solver_data["D"]
            raise InvalidSystemMatrixException("Invalid diffusion matrix. The sum of the columns does not sum to zero. " + str(maxcolsum) + ' ' + str(colsum[0,maxcolsum]) + "\nThis can be caused by a large difference between the largest and smallest diffusion coefficients.")

    def create_connectivity_matrix(self):
        """ Assemble a connectivity matrix in CCS format. """

        fs = self.mesh.get_function_space()
        trial_function = dolfin.TrialFunction(fs)
        test_function = dolfin.TestFunction(fs)
        a_K = -1*dolfin.inner(dolfin.nabla_grad(trial_function), dolfin.nabla_grad(test_function)) * dolfin.dx
        K = dolfin.assemble(a_K)
        rows, cols, vals = dolfin.as_backend_type(K).data()
        Kcrs = scipy.sparse.csc_matrix((vals, cols, rows))
        return Kcrs

    def get_solver_datastructure(self):
        """ Return the datastructures needed by the URDME solvers.
           get_solver_datastructure() creates and populates a dictionary, urdme_solver_data,
           containing the mandatory input data structures of the core NSM solver in URDME.
           Those data structures are
           D    - the Diffusion matrix
           N    - the stochiometry matrix
           G    - the dependency graph
           vol  - the volume vector
           sd   - the subdomain vector
           data - the data vector
           u0   - the intial condition
           The follwing data is also returned, unlike in the legacy URDME interface:
           p - the vertex coordinates
           K - a (Nvoxel x Nvoxel) connectivity matrix
        """
        
        urdme_solver_data = {}
        num_species = self.get_num_species()

        # Stoichimetric matrix
        N = self.create_stoichiometric_matrix()
        urdme_solver_data['N'] = N
        # Dependency Graph
        G = self.create_dependency_graph()
        urdme_solver_data['G']  = G

        # Volume vector
        result =  self.create_system_matrix()
        vol = result['vol']
        urdme_solver_data['dofvolumes'] = vol

        #TODO: Make use of all dofs values, requires modification of the core URDME solver.
        self.dofvol = vol[::len(self.listOfSpecies)]
        urdme_solver_data['vol'] = self.dofvol

        D = result['D']
        urdme_solver_data['D'] = D

        num_dofvox = self.dofvol.shape[0]

        # Get vertex to dof ordering
        vertex_to_dof = dolfin.vertex_to_dof_map(self.mesh.get_function_space())
        dof_to_vertex = dolfin.dof_to_vertex_map(self.mesh.get_function_space())

        vertex_to_dof_to_vertex = dof_to_vertex[vertex_to_dof]

        # Subdomain vector
        # convert to dof ordering
        sd_vec_dof = numpy.zeros(num_dofvox)
        for ndx, sd_val in enumerate(self.get_subdomain_vector()):
            sd_vec_dof[vertex_to_dof[ndx]] = sd_val
        urdme_solver_data['sd'] = sd_vec_dof

        # Data vector. If not present in model, it defaults to a vector with all elements zero.
        # convert to dof ordering
        data = numpy.zeros((1, num_dofvox))
        if len(self.listOfDataFunctions) > 0:
            data = numpy.zeros((len(self.listOfDataFunctions), num_dofvox))
            coords = self.mesh.coordinates()
            for ndf, df in enumerate(self.listOfDataFunctions):
                for ndx in range(len(coords)):
                    vox_coords = numpy.zeros(3)
                    for cndx in range(len(coords[ndx])):
                        vox_coords[cndx] = coords[ndx][cndx]
                    data[ndf][vertex_to_dof[ndx]] = df.map(vox_coords)

        urdme_solver_data['data'] = data

        if not hasattr(self,'u0'):
            self.initialize_initial_condition()

        # Initial Conditions, convert to dof ordering
        u0_dof = numpy.zeros((num_species, num_dofvox))
        for vox_ndx in range(self.mesh.get_num_voxels()):
            dof_ndx = vertex_to_dof[vox_ndx]
            # With periodic BCs the same dof_ndx voxel will get written to twice
            # which may overwrite the value.  We need to check for this case.
            if vertex_to_dof_to_vertex[vox_ndx] != vox_ndx:
                vox_ndx2 = vertex_to_dof_to_vertex[vox_ndx]
                for cndx in range(num_species):
                    if self.u0[cndx, vox_ndx] == 0 or self.u0[cndx, vox_ndx] == self.u0[cndx, vox_ndx2]:
                        u0_dof[cndx, dof_ndx] = self.u0[cndx, vox_ndx]
                    elif self.u0[cndx, vox_ndx2] == 0 and vox_ndx < vox_ndx2:
                        self.u0[cndx, vox_ndx2] = self.u0[cndx, vox_ndx]
                    elif self.u0[cndx, vox_ndx2] == 0 and vox_ndx > vox_ndx2:
                        u0_dof[cndx, dof_ndx] = self.u0[cndx, vox_ndx]
                    else:
                        sys.stderr.write("Warning: the initial condition for species {0} in voxel {1} will be discarded due to periodic boundary conditions.\n".format(self.listOfSpecies.keys()[cndx], vox_ndx))
            else:
                for cndx in range(num_species):
                    u0_dof[cndx, dof_ndx] = self.u0[cndx, vox_ndx]
        urdme_solver_data['u0'] = u0_dof

        tspan = numpy.asarray(self.tspan, dtype=numpy.float)
        urdme_solver_data['tspan'] = tspan

        # Vertex coordinates
        # convert to dof ordering
        p_dof = numpy.zeros((num_dofvox, 3))
        for vox_ndx, row in enumerate(self.mesh.get_voxels()):
            p_dof[vertex_to_dof[vox_ndx], :len(row)] = row
        urdme_solver_data['p'] = p_dof

        # Connectivity matrix
        urdme_solver_data['K'] = self.create_connectivity_matrix()

        urdme_solver_data['report']=0

        return urdme_solver_data


    def assemble(self):
        """  Assemble the mass and stiffness matrices using Dolfin.
            Returns: A dictionary containing two dictionaries, one for the stiffness matrices
            and one for the mass matrices. Those dictionaries has the species names as keys and
            the matrices are in CSR format.
            """

        if self.xmesh == None:
            self.create_extended_mesh()

        self._initialize_species_to_subdomains()

        function_space = self.xmesh.function_space
        trial_functions = OrderedDict()
        test_functions = OrderedDict()
        stiffness_matrices = OrderedDict()
        mass_matrices = OrderedDict()

        # The maximum dimension that a species is active on (currently not used)
        maxdim = 1
        for spec in self.listOfSpecies:
            dim = self.listOfSpecies[spec].dim()
            if dim > maxdim:
                maxdim = dim

        for spec in self.listOfSpecies:
            trial_functions[spec] = dolfin.TrialFunction(function_space[spec])
            test_functions[spec] = dolfin.TestFunction(function_space[spec])


        weak_form_K = {}
        weak_form_M = {}

        ndofs = None

        # Set up the forms
        for spec_name, species in self.listOfSpecies.items():

            # Find out what subdomains this species is active on
            weak_form_K[spec_name] = dolfin.inner(dolfin.nabla_grad(trial_functions[spec_name]), dolfin.nabla_grad(test_functions[spec_name]))*dolfin.dx
            weak_form_M[spec_name] = trial_functions[spec_name]*test_functions[spec_name]*dolfin.dx

        # Assemble the matrices
        for spec_name, species in self.listOfSpecies.items():
            stiffness_matrices[spec_name] = dolfin.assemble(weak_form_K[spec_name])
            if ndofs is None:
                ndofs = stiffness_matrices[spec_name].size(0)
                self.mesh.set_num_dof_voxels(ndofs)

            # We cannot include the diffusion constant in the assembly, dolfin does not seem to deal well
            # with small diffusion constants (it drops small elements), so we multiply here. 
            stiffness_matrices[spec_name] = species.diffusion_constant * stiffness_matrices[spec_name]
            mass_matrices[spec_name] = dolfin.assemble(weak_form_M[spec_name])


        return {'K':stiffness_matrices, 'M':mass_matrices}

    def run(self, number_of_trajectories=1, solver='nsm', seed=None, report_level=0):
        """ Simulate the model.
        Args:
            solver: A str or class type that is a subclass of URDMESolver.  Default: NSM solver.
            number_of_trajectories: How many trajectories should be run.
            seed: An int, the random seed given to the solver.
            report_level: An int, Level of output from the solver: 0, 1, or 2. Default: 0.
        Returns:
            A URDMEResult object with the results of the simulation.
        """

        #If solver is a subclass of URDMESolver, use it directly.
        if isinstance(solver, (type, types.ClassType)) and  issubclass(solver, URDMESolver):
            sol = solver(self, report_level=report_level)
        elif type(solver) is str:
            if solver == 'nsm':
                from nsmsolver import NSMSolver
                sol = NSMSolver(self, report_level=report_level)
            else:
                raise URDMEError("Unknown solver: {0}".format(solver_name))
        else:
            raise URDMEError("solver argument to urdme() must be a string or a URDMESolver class object.")

        return sol.run(number_of_trajectories=number_of_trajectories, seed=seed)


class URDMEMesh(dolfin.Mesh):
    """ A URDME mesh extends the Dolfin mesh class.
        Provides wrappers around dolfins built-in simple geometries/mesh generation function.
        These following methods will all give regular meshes that will produce discretizations that are equivalent to Cartesian grids.
    """

    def __init__(self, mesh=None):
        self.constrained_domain = None
        dolfin.Mesh.__init__(self, mesh)
        self.function_space = None
        self.num_dof_voxels = None
        self.init()


    def __getstate__(self):

        state = {}
        state['function_space'] = None

        tmpfile = tempfile.NamedTemporaryFile(suffix=".xml")
        dolfin.File(tmpfile.name) << self
        tmpfile.seek(0)

        state['meshdata'] = tmpfile.read()
        tmpfile.close()

        if self.constrained_domain is not None:
            # Warning: This is black magic.
            try:
                cdd = {}
                cdd['source'] = inspect.getsource(self.constrained_domain.__class__)
                cdd['name'] = self.constrained_domain.__class__.__name__
                cdd['dict'] = {}
                for k, v in self.constrained_domain.__dict__.iteritems():
                    if type(v).__name__ != 'SwigPyObject':
                        cdd['dict'][k] = v
                state['constrained_domain'] = cdd
            except Exception as e:
                sys.stderr.write("error pickling mesh.constrained_domain: {0}\n".format(e))
                raise e
        if self.num_dof_voxels is not None:
            state['num_dof_voxels'] = self.num_dof_voxels

        return state

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """

        try:
            fd = tempfile.NamedTemporaryFile(suffix=".xml")
            fdname = fd.name
            fd.write(state['meshdata'])
            fd.seek(0)
            self.__init__(fd.name)

            if 'constrained_domain' in state and state['constrained_domain'] is not None:
                # Black magic to match that in __getstate__
                cdd = state['constrained_domain']
                compiled_class = compile(cdd['source'], 'pyurdme.mesh.constrained_domain', 'exec')
                eval(compiled_class)
                compiled_object = eval("{0}()".format(cdd['name']))
                for k, v in cdd['dict'].iteritems():
                    compiled_object.__dict__[k] = v
                self.constrained_domain = compiled_object
            if 'num_dof_voxels' in state and state['num_dof_voxels'] is not None:
                self.num_dof_voxels = state['num_dof_voxels']
        except Exception as e:
            print "Error unpickling model, could not recreate the mesh."
            raise e

    def add_periodic_boundary_condition(self, domain):
        """ Add a periodic boundary mapping object (a subclass of dolfin.SubDomain). """
        self.constrained_domain = domain

    def get_function_space(self):
        """ Get the FunctionSpace dolfin object for this mesh. """
        if self.function_space is not None:
            return self.function_space
        else:
            if self.constrained_domain is not None:
                fs = dolfin.FunctionSpace(self, "Lagrange", 1, constrained_domain=self.constrained_domain)
            else:
                fs = dolfin.FunctionSpace(self, "Lagrange", 1)
            self.function_space = fs
            return fs

    def get_num_voxels(self):
        """ Get the number of voxels in the vertex ordering. """
        return self.num_vertices()

    def set_num_dof_voxels(self, num):
        """ Set the number of voxels in the DOF ordering. """
        self.num_dof_voxels = num

    def get_num_dof_voxels(self):
        """ Get the number of voxels in the DOF ordering. """
        if self.num_dof_voxels is None:
            raise URDMEError('NumDofVoxels is not set')
        return self.num_dof_voxels

    def get_voxels(self):
        """ Return the (x,y,z) coordinate of each voxel. """
        coords = self.coordinates()
        if coords.shape[1] == 2:
            coords = numpy.append(coords, numpy.tile([0],(coords.shape[0],1)), 1)
        return coords

    def closest_vertex(self, x):
        """ Get index of the vertex in the coordinate list closest to the point x. """
        coords = self.get_voxels()
        shape = coords.shape

        if isinstance(x, (int, float)):
            x = [x]

        if len(x) == 2:
            point = numpy.append(x, 0.0)
        else:
            point = x
        reppoint = numpy.tile(point, (shape[0], 1))
        dist = numpy.sqrt(numpy.sum((coords-reppoint)**2, axis=1))
        ix = numpy.argmin(dist)
        return ix

    def get_mesh_size(self):
        """ Estimate of mesh size at each vertex. """
        coordinates = self.coordinates()

        # Compute the circumradius of the cells
        cr = []
        for i in range(self.num_cells()):
            cell = dolfin.Cell(self, i)
            cr.append(cell.diameter()/2.0)

        # Compute the mean for each vertex based on all incident cells
        vtx2cell = self.topology()(0,self.topology().dim())
        vtxh = []
        for i in range(self.num_vertices()):
            v2c = vtx2cell(i)
            h = 0.0
            for indx in v2c:
                h += cr[indx]
            h = h/len(v2c)
            vtxh.append(h)

        return vtxh

    def get_normalized_coordinates(self):
        """ Return vertex coordinates centered at origin. """

        # Compute mesh centroid
        vtx = self.coordinates()
        centroid = numpy.mean(vtx,axis=0)
        # Shift so the centroid is now origo
        normalized_vtx = numpy.zeros(numpy.shape(vtx))
        for i,v in enumerate(vtx):
            normalized_vtx[i,:] = v - centroid

        return normalized_vtx

    def get_scaled_normalized_coordinates(self):
        """ Return vertex coordinates scaled to the interval (-1,1) and centered at origin. """
        # Scale the verices so the max dimension is in the range (-1,1) to be compatible with the browser display
        vtx = self.coordinates()
        maxvtx = numpy.max(numpy.amax(vtx,axis=0))
        factor = 1/maxvtx
        vtx = factor*vtx

        # Compute mesh centroid
        centroid = numpy.mean(vtx,axis=0)
        # Shift so the centroid is now origo
        normalized_vtx = numpy.zeros(numpy.shape(vtx))
        for i,v in enumerate(vtx):
            normalized_vtx[i,:] = v - centroid


        return factor, normalized_vtx

    def get_scaled_coordinates(self):
        """ Return vertex coordinates scaled to the interval (-1,1). """
        # Scale the verices so the max dimension is in the range (-1,1) to be compatible with the browser display
        vtx = self.coordinates()
        maxvtx = numpy.max(numpy.amax(vtx,axis=0))
        factor = 1/maxvtx
        return factor, factor*vtx

    @classmethod
    def fromPoints(cls, points):
        """ Create a mesh from a list of points (3D) only. Points is a list or numpy array
            
            [[x1,y1,z1],
             [x2,y2,z2],
              ...
              ...]
        
        """
        try:
            import lxml.etree as etree
            no_pretty_print = False
        except:
            import xml.etree.ElementTree as etree
            import xml.dom.minidom
            import re
            no_pretty_print = True
        
        # Create a Delauny triangulation of the points
        import scipy.spatial
        tri = scipy.spatial.Delaunay(points, furthest_site=False)
        #tri = scipy.spatial.Delaunay(points)
        
        # Write a temporary Dolfin XML file.
        tree = etree.Element('dolfin')
        mesh = etree.Element('mesh')
        mesh.set('celltype', 'tetrahedron')
        mesh.set('dim', '3')
        vertices = etree.Element('vertices')
        dim = numpy.shape(tri.points)
        vertices.set('size',str(dim[0]))

        for i,v in enumerate(tri.points):
            vtx = etree.Element('vertex')
            vtx.set('index',str(i))
            vtx.set('x',str(v[0]))
            vtx.set('y',str(v[1]))
            vtx.set('z',str(v[2]))
            vertices.append(vtx)

        mesh.append(vertices)
        dim = numpy.shape(tri.simplices)
        cells = etree.Element('cells')
        cells.set('size',str(dim[0]))
        for i,cell in enumerate(tri.simplices):
            c = etree.Element('tetrahedron')
            c.set('index',str(i))
            c.set('v0',str(cell[0]))
            c.set('v1',str(cell[1]))
            c.set('v2',str(cell[2]))
            c.set('v3',str(cell[3]))
            cells.append(c)
        mesh.append(cells)
        tree.append(mesh)

        f = tempfile.NamedTemporaryFile(suffix='.xml', delete=False)
        filename = f.name
        with open(f.name,'w') as fh:
            fh.write(etree.tostring(tree))
        msh = URDMEMesh(dolfin.Mesh(filename))
        os.remove(filename)
        return msh
            

    @classmethod
    def generate_unit_interval_mesh(cls, nx, periodic=False):
        """ Unit Interval (1D) of with nx points in the axes. """
        return cls.generate_interval_mesh(nx=nx, a=0, b=1, periodic=periodic)

    @classmethod
    def generate_unit_square_mesh(cls, nx, ny, periodic=False):
        """ Unit Square (2D) of with nx, ny points in the respective axes. """
        return cls.generate_square_mesh(L=1, nx=nx, ny=ny, periodic=periodic)

    @classmethod
    def generate_unit_cube_mesh(cls, nx, ny, nz, periodic=False):
        """ Unit Cube (3D) of with nx, ny, nz points in the respective axes. """
        return cls.generate_cube_mesh(L=1,nx=nx, ny=ny, nz=nz, periodic=periodic)

    @classmethod
    def generate_interval_mesh(cls, nx, a, b, periodic=False):
        """ Interval (1D) of with nx points in the axes, and side length L. """
        mesh = dolfin.IntervalMesh(nx, a, b)
        ret = URDMEMesh(mesh)
        if isinstance(periodic, bool) and periodic:
            ret.add_periodic_boundary_condition(IntervalMeshPeriodicBoundary(a=a, b=b))
        elif isinstance(periodic, dolfin.SubDomain):
            ret.add_periodic_boundary_condition(periodic)
        return ret

    @classmethod
    def generate_square_mesh(cls, L, nx, ny, periodic=False):
        """ Unit Square (2D) of with nx, ny points in the respective axes, and side length L. """
        try:
            mesh = dolfin.RectangleMesh(0, 0, L, L, nx, ny)
        except (TypeError, NotImplementedError) as e:
            # for Dolfin 1.6+
            rect = mshr.Rectangle(dolfin.Point(0,0), dolfin.Point(L,L))
            mesh = mshr.generate_mesh(rect, nx)
        ret = URDMEMesh(mesh)
        if isinstance(periodic, bool) and periodic:
            ret.add_periodic_boundary_condition(SquareMeshPeriodicBoundary(Lx=L, Ly=L))
        elif isinstance(periodic, dolfin.SubDomain):
            ret.add_periodic_boundary_condition(periodic)
        return ret

    @classmethod
    def generate_cube_mesh(cls, L, nx, ny, nz, periodic=False):
        """ Unit Cube (3D) of with nx, ny, nz points in the respective axes, and side length L. """
        try:
            mesh = dolfin.BoxMesh(0, 0, 0, L, L, L, nx, ny, nz)
        except (TypeError, NotImplementedError) as e:
            # for Dolfin 1.6+
            box = mshr.Box(dolfin.Point(0,0,0), dolfin.Point(L,L,L))
            mesh = mshr.generate_mesh(box, nx)
        ret = URDMEMesh(mesh)
        if isinstance(periodic, bool) and periodic:
            ret.add_periodic_boundary_condition(CubeMeshPeriodicBoundary(Lx=L, Ly=L, Lz=L))
        elif isinstance(periodic, dolfin.SubDomain):
            ret.add_periodic_boundary_condition(periodic)
        return ret

    @classmethod
    def read_dolfin_mesh(cls, filename=None, colors = []):
        """ Import a mesh in Dolfins native .xml format """

        try:
            dolfin_mesh = dolfin.Mesh(filename)
            mesh = URDMEMesh(mesh=dolfin_mesh)
            return mesh
        except Exception as e:
            raise MeshImportError("Failed to import mesh: " + filename+"\n" + str(e))

    @classmethod
    def read_mesh(cls, filename=None, colors = []):
        """Import a mesh in gmsh .msh or Dolfins .xml format"""
        if filename[-4:]==".msh":
            #if the input file is a .msh, we convert it into a Dolfin .xml
            subprocess.call(["dolfin-convert",filename,filename[:-4]+".xml"])
        mesh = cls.read_dolfin_mesh(filename[:-4]+".xml",colors)
        return mesh

    @classmethod
    def read_geometry(cls, filename=None, dimension=2, clscale=1, colors=[]):
        """Import a mesh from a geometry"""
        mesh_filename = (filename[:-4] if filename[-4:]==".geo" else filename)+".msh"
        subprocess.call(["gmsh","-"+str(dimension),"-clscale",str(clscale),filename,"-o",mesh_filename])
        mesh = cls.read_mesh(mesh_filename,colors)
        return mesh

    def export_to_three_js(self, colors = None):
        """ return a Json string of the mesh in THREE Js format.
            If a colors list is specified, it should have the num_voxels entries
        """
        self.init(2,0)
        document = {}
        document["metadata"] = {"formatVersion":3}
        gfdg,vtx = self.get_scaled_normalized_coordinates()

        if self.topology().dim() == 2:
            # 2D
            num_elements = self.num_cells()
            # This is a fix for the built-in 2D meshes that only have x,y-coordinates.
            dims = numpy.shape(vtx)
            if dims[1] == 2:
                vtxx = numpy.zeros((dims[0],3))
                for i, v in enumerate(vtx):
                    vtxx[i,:]=(list(v)+[0])
                vtx = vtxx
        else:
            # 3D
            num_elements = self.num_facets()

        materials = [ {
                     "DbgColor" : 15658734,
                     "DbgIndex" : 0,
                     "DbgName" : "dummy",
                     "colorDiffuse" : [ 1, 0, 0 ],
                     } ]

        document["materials"] = materials
        document["vertices"] = list(vtx.flatten())

        if colors == None:
            # Default color is blue
            colors = [255]*self.num_vertices()

        document["colors"] = colors

        connectivity = self.topology()(2,0)
        faces = []

        for i in range(num_elements):
            face = connectivity(i)
            f = []
            for ind in face:
                if int(ind) >= self.num_vertices():
                    raise Exception("Out of bounds")

                f.append(int(ind))
            faces += ([128]+f+f)
        document["faces"] = list(faces)

        #Test that we can index into vertices
        vertices = document["vertices"]

        return json.dumps(document)

    def _ipython_display_(self, filename=None, colors=None, width=500):
        self.display(filename=filename, colors=colors, width=width)

    def display(self, filename=None, colors=None, width=500, camera=[0,0,1]):
        load_pyurdme_javascript_libraries()
        jstr = self.export_to_three_js(colors=colors)
        hstr = None
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/mesh.html",'r') as fd:
            hstr = fd.read()
        if hstr is None:
            raise Exception("could note open template mesh.html")
        hstr = hstr.replace('###PYURDME_MESH_JSON###',jstr)
        hstr = hstr.replace('###WIDTH###',str(width))
        height = int(width * 0.75)
        # Create a random id for the display div. This is to avioid multiple plots ending up in the same
        # div in Ipython notebook
        displayareaid=str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###',displayareaid)
        # ###CAMERA_X###, ###CAMERA_Y###, ###CAMERA_Z###
        hstr = hstr.replace('###CAMERA_X###',str(camera[0]))
        hstr = hstr.replace('###CAMERA_Y###',str(camera[1]))
        hstr = hstr.replace('###CAMERA_Z###',str(camera[2]))
        html = '<div style="width: {0}px; height: {1}px;" id="{2}" ></div>'.format(width, height, displayareaid)

        if filename is not None:
            with open(filename, 'w') as fd:
                fd.write("""
<html>
    <head>
        <title>PyURDME Result</title> <style>canvas { width: 100%; height: 100% }</style> </head>
        <body>
""")
                fd.write(html+hstr)
                fd.write("""
        </body>
</html>""")
        else:
            IPython.display.display(IPython.display.HTML(html+hstr))



class URDMEXmesh():
    """ Extended mesh object.
        Contains function spaces and dof mappings.
    """

    def __init__(self):
        self.coordinates = None
        self.function_space = {}
        self.vertex_to_dof_map = {}
        self.dof_to_vertex_map = {}

