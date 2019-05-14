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
    raise Exception("SpatialPy requires h5py.")

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
#dolfin.set_log_level(dolfin.ERROR)
#import logging
#logging.getLogger('FFC').setLevel(logging.ERROR)
#logging.getLogger('UFL').setLevel(logging.ERROR)

class SPATIALPYModel(Model):
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


