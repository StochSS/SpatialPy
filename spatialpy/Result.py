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

from spatialpy.VTKReader import VTKReader
from spatialpy.Model import *

import inspect

import pickle
import json


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



class Result(dict):
    """ Result object for a URDME simulation, extends the dict object. """

    def __init__(self, model=None, result_dir=None, loaddata=False):
        self.model = model
        self.U = None
        self.tspan = None
        self.data_is_loaded = False
        self.stdout = None
        self.stderr = None
        self.result_dir = result_dir



#    def get_endtime_model(self):
#        """ Return a URDME model object with the initial conditions set to the final time point of the
#            result object.
#        """
#        if self.model is None:
#            raise Exception("can not continue a result with no model")
#        # create a soft copy
#        model_str = pickle.dumps(self.model)
#        model2 = pickle.loads(model_str)
#        # set the initial conditions
#        model2.u0 = numpy.zeros(self.model.u0.shape)
#        for s, sname in enumerate(self.model.listOfSpecies):
#            model2.u0[s,:] = self.get_species(sname, timepoints=-1)
#        return model2


#    def __getstate__(self):
#        """ Used by pickle to get state when pickling. We need to read the contents of the
#        output file since we can't pickle file objects. """
#
#        try:
#            with open(self.filename,mode='rb') as fh:
#                filecontents = fh.read()
#        except Exception as e:
#            raise Exception(("Error pickling model. Failed to read result file:",str(e)))
#
#        state = self.__dict__
#        state["filecontents"] = filecontents
#
#        state["v2d"] = self.get_v2d()
#        state["d2v"] = self.get_d2v()
#
#        return state
#
#
#    def __setstate__(self, state):
#        """ Used by pickle to set state when unpickling. """
#
#        # If the object contains filecontents, write those to a new tmp file.
#        try:
#            filecontents = state.pop("filecontents",None)
#            fd = tempfile.NamedTemporaryFile(delete=False, dir=os.environ.get('PYURDME_TMPDIR'))
#            with open(fd.name, mode='wb') as fh:
#                fh.write(filecontents)
#            state["filename"] = fd.name
#        except Exception as e:
#            print "Error unpickling model, could not recreate the solution file."
#            raise e
#
#        for k,v in state.items():
#            self.__dict__[k] = v


    def read_step(self, step_num):
        """ Read the data for simulation step 'step_num'. """
        reader = VTKReader()
        filename = os.path.join(self.result_dir, "output{0}.vtk".format(step_num))
        #print("read_step({0}) opening '{1}'".format(step_num, filename))
        reader.setfilename(filename)
        reader.readfile()
        if reader.getpoints() is None or reader.getarrays() is None:
            raise ResultError("read_step(step_num={0}): got data = None".format(step_num))
        points = reader.getpoints()
        vtk_data = reader.getarrays()
        return (points, vtk_data)


#    def read_solution(self):
#        """ Read the tspan and U matrix into memory. """
#        raise Exception("todo")
##        self.U = U
##        self.tspan = tspan
##        self.data_is_loaded = True

    def get_timespan(self):
        self.tspan = numpy.linspace(0,self.model.num_timesteps,
                num=self.model.num_timesteps+1) * self.model.timestep_size
        return self.tspan

    def get_species(self, species, timepoints=None, concentration=False):
        """ Get the populations/concentration values for a given species in the model for 
            one or all timepoints.  
            
            If 'timepoints' is None (default), a matrix of dimension:
            (number of timepoints) x (number of voxels) is returned.  If an integer value is
            given, that value is used to index into the timespan, and that time point is returned
            as a 1D array with size (number of voxel). 

            If concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.
        """


        species_map = self.model.species_map
        num_species = self.model.get_num_species()
        num_voxel = self.model.mesh.get_num_voxels()

        if isinstance(species,str):
            spec_name = species
        else:
            spec_name = species.name

        if spec_name not in self.model.listOfSpecies.keys():
            raise ResultError("Species '{0}' not found".format(spec_name))

        t_index_arr = numpy.linspace(0,self.model.num_timesteps, 
                            num=self.model.num_timesteps+1, dtype=int)

        if timepoints is not None:
            if isinstance(timepoints,float):
                raise ResultError("timepoints argument must be an integer, the index of time timespan")
            t_index_arr = t_index_arr[timepoints]
        #if not isinstance(t_index_arr,list):
        #    t_index_arr = [t_index_arr]
        try:
            num_timepoints = len(t_index_arr)
        except Exception as e:
            t_index_arr = [t_index_arr]
            num_timepoints = 1

        ret = numpy.zeros( (num_timepoints, num_voxel))
        for ndx, t_ndx in enumerate(t_index_arr):
            (_, step) = self.read_step(t_ndx)
            if concentration:
                # concentration = (copy_number/volume)
                # volume = (mass/density)
                ret[ndx,:] = step['C['+spec_name+']'] / (step['mass'] / step['rho'] )
            else:
                ret[ndx,:] = step['C['+spec_name+']']
        if ret.shape[0] == 1:
            ret = ret.flatten()
        return ret


#    def __setattr__(self, k, v):
#        if k in self.keys():
#            self[k] = v
#        elif not hasattr(self, k):
#            self[k] = v
#        else:
#            raise AttributeError("Cannot set '%s', cls attribute already exists" % ( k, ))
#
#    def __setupitems__(self, k):
#        if (k == 'U' or k == 'tspan') and not self.data_is_loaded:
#            if self.result_dir is None:
#                raise AttributeError("This result object has no data file.")
#            self.read_solution()
#
#    def __getitem__(self, k):
#        self.__setupitems__(k)
#        if k in self.keys():
#            return self.get(k)
#        raise KeyError("Object has no attribute {0}".format(k))
#
#    def __getattr__(self, k):
#        self.__setupitems__(k)
#        if k in self.keys():
#            return self.get(k)
#        raise AttributeError("Object has no attribute {0}".format(k))

    def __del__(self):
        """ Deconstructor. """
        #   if not self.data_is_loaded:
        try:
            if self.result_dir is not None:
                try:
                    shutil.rmtree(self.result_dir)
                except OSError as e:
                    print("Could not delete '{0}'".format(self.result_dir))
        except Exception as e:
            pass




    def export_to_csv(self, folder_name):
        """ Dump trajectory to a set CSV files, the first specifies the mesh (mesh.csv) and the rest specify trajectory data for each species (species_S.csv for species named 'S').
            The columns of mesh.csv are: 'Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Subdomain'.
            The columns of species_S.csv are: 'Time', 'Voxel 0', Voxel 1', ... 'Voxel N'.
        """
        import csv
        subprocess.call(["mkdir", "-p", folder_name])
        #['Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Subdomain']
        with open(os.path.join(folder_name,'mesh.csv'), 'w+') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Subdomain'])
            vol = self.model.get_solver_datastructure()['vol']
            for ndx in range(self.model.mesh.get_num_voxels()):
                row = [ndx]+self.model.mesh.coordinates()[ndx,:].tolist()+[vol[ndx]]+[self.model.mesh.sd[ndx]]
                writer.writerow(row)

        for spec in self.model.listOfSpecies:
            #['Time', 'Voxel 0', Voxel 1', ... 'Voxel N']
            with open(os.path.join(folder_name,'species_{0}.csv'.format(spec)), 'w+') as csvfile:
                data = self.get_species(spec)
                (num_t,num_vox) = data.shape
                writer = csv.writer(csvfile, delimiter=',')
                row = ['Time']
                for v in range(num_vox):
                    row.append('Voxel {0}'.format(v))
                writer.writerow(row)
                timespan = self.get_timespan()
                for t in range(num_t):
                    writer.writerow([timespan[t].tolist()] + data[t,:].tolist())

    def export_to_vtk(self, species, folder_name):
        """ Dump the trajectory to a collection of vtk files in the folder folder_name (created if non-existant).
            The exported data is #molecules/volume, where the volume unit is implicit from the mesh dimension. """

#        #self._initialize_sol()
#        subprocess.call(["mkdir", "-p", folder_name])
#        fd = dolfin.File(os.path.join(folder_name, "trajectory.xdmf").encode('ascii', 'ignore'))
#        func = dolfin.Function(self.model.mesh.get_function_space())
#        func_vector = func.vector()
#        vertex_to_dof_map = self.get_v2d()
#
#        for i, time in enumerate(self.tspan):
#            solvector = self.get_species(species,i,concentration=True)
#            for j, val in enumerate(solvector):
#                func_vector[vertex_to_dof_map[j]] = val
#            fd << func
        raise Exception("todo")

    def export_to_xyx(self, filename, species=None, file_format="VMD"):
        """ Dump the solution attached to a model as a xyz file. This format can be
            read by e.g. VMD, Jmol and Paraview. """

        if self.U is None:
            raise URDMEError("No solution found in the model.")

        #outfile = open(filename,"w")
        dims = numpy.shape(self.U)
        Ndofs = dims[0]
        Mspecies = len(self.model.listOfSpecies)
        Ncells = Ndofs / Mspecies

        coordinates = self.model.mesh.get_voxels()
        coordinatestr = coordinates.astype(str)

        if species == None:
            species = list(self.model.listOfSpecies.keys())

        if file_format == "VMD":
            outfile = open(filename, "w")
            filestr = ""
            for i, time in enumerate(self.tspan):
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(species):
                    for k in range(Ncells):
                        for mol in range(self.U[k * Mspecies + j, i]):
                            # Sample a random position in a sphere of radius computed from the voxel volume
                            # TODO: Sample volume
                            linestr = spec + "\t" + '\t'.join(coordinatestr[k, :]) + "\n"
                            filestr += linestr

            outfile.write(filestr)
            outfile.close()

        elif file_format == "ParaView":
            foldername = filename
            os.mkdir(foldername)
            for i, time in enumerate(self.tspan):
                outfile = open(foldername + "/" + filename + "." + str(i), "w")
                number_of_atoms = numpy.sum(self.U[:, i])
                filestr = ""
                filestr += (str(number_of_atoms) + "\n" + "timestep " + str(i) + " time " + str(time) + "\n")
                for j, spec in enumerate(self.model.listOfSpecies):
                    for k in range(Ncells):
                        for mol in range(model.U[k * Mspecies + j, i]):
                            linestr = spec + "\t" + '\t'.join(coordinatestr[k, :]) + "\n"
                            filestr += linestr
                outfile.write(filestr)
                outfile.close()



    def _export_to_particle_js(self,species,time_index, colors=None):
        """ Create a html string for displaying the particles as small spheres. """
        import random
        with open(os.path.dirname(os.path.abspath(__file__))+"/data/three.js_templates/particles.html",'r') as fd:
            template = fd.read()

        factor, coordinates = self.model.mesh.get_scaled_normalized_coordinates()
        dims = numpy.shape(coordinates)
        if dims[1]==2:
            is3d = 0
            vtxx = numpy.zeros((dims[0],3))
            for i, v in enumerate(coordinates):
                vtxx[i,:]=(list(v)+[0])
            coordinates = vtxx
        else:
            is3d = 1

        h = self.model.mesh.get_mesh_size()

        x=[]
        y=[]
        z=[]
        c=[]
        radius = []

        total_num_particles = 0
        #colors = ["blue","red","yellow", "green"]
        if colors == None:
            colors =  get_N_HexCol(len(species))

        if not isinstance(species, list):
           species = [species]
        for j,spec in enumerate(species):
            timeslice = self.get_species(spec, time_index)
            ns = numpy.sum(timeslice)
            total_num_particles += ns

            for i, particles in enumerate(timeslice):
                # "Radius" of voxel
                hix = h[i]*factor
                hiy = hix
                hiz = hix*is3d

                for particle in range(int(particles)):
                    x.append((coordinates[i,0]+random.uniform(-1,1)*hix))
                    y.append((coordinates[i,1]+random.uniform(-1,1)*hiy))
                    z.append((coordinates[i,2]+random.uniform(-1,1)*hiz))
                    if self.model.listOfSpecies[spec].reaction_radius:
                        radius.append(factor*self.model.listOfSpecies[spec].reaction_radius)
                    else:
                        radius.append(0.01)

                    c.append(colors[j])

        template = template.replace("__X__",str(x))
        template = template.replace("__Y__",str(y))
        template = template.replace("__Z__",str(z))
        template = template.replace("__COLOR__",str(c))
        template = template.replace("__RADIUS__",str(radius))

        return template


    def export_to_three_js(self, species, time_index):
        """ Return a json serialized document that can
            be read and visualized by three.js.
        """

        colors = self._compute_solution_colors(species,time_index)
        return self.model.mesh.export_to_three_js(colors=colors)

    def _copynumber_to_concentration(self,copy_number_data):
        """ Scale comy numbers to concentrations (in unit mol/volume),
            where the volume unit is defined by the user input.
            Dof-ordering is assumed in both solution and volumes.
        """

        shape = numpy.shape(copy_number_data)
        if len(shape) == 1:
            shape = (1,shape[0])

        scaled_sol = numpy.zeros(shape)
        scaled_sol[:,:] = copy_number_data
        dims = numpy.shape(scaled_sol)

        for t in range(dims[0]):
            timeslice = scaled_sol[t,:]
            for i,cn in enumerate(timeslice):
                scaled_sol[t, i] = float(cn)/(6.022e23*self.model.vol[i])

        return scaled_sol

    @classmethod
    def _compute_colors(cls, x):
        import matplotlib.cm

        # Get RGB color map proportional to the concentration.
        cm = matplotlib.cm.ScalarMappable()
        crgba= cm.to_rgba(x, bytes = True)
        # Convert RGB to HEX
        colors= []
        for row in crgba:
            # get R,G,B of RGBA
            colors.append(_rgb_to_hex(tuple(list(row[0:3]))))

        # Convert Hex to Decimal
        for i,c in enumerate(colors):
            colors[i] = int(c,0)


        return colors


    def _compute_solution_colors(self,species, time_index):
        """ Create a color list for species at time. """

        timeslice = self.get_species(species,time_index, concentration = True)
        colors = self._compute_colors(timeslice)
        return colors

    def display_particles(self,species, time_index, width=500):
        load_pyurdme_javascript_libraries()
        hstr = self._export_to_particle_js(species, time_index)
        displayareaid=str(uuid.uuid4())
        hstr = hstr.replace('###DISPLAYAREAID###',displayareaid)
        hstr = hstr.replace('###WIDTH###',str(width))
        height = int(width*0.75)

        html = '<div style="width: {0}px; height: {1}px;" id="{2}" ></div>'.format(width, height, displayareaid)
        IPython.display.display(IPython.display.HTML(html+hstr))


    def display(self, species, time_index, opacity=1.0, wireframe=True, width=500, camera=[0,0,1]):
        """ Plot the trajectory as a PDE style plot. """
#        load_pyurdme_javascript_libraries()
#        data = self.get_species(species,time_index,concentration=True)
#        fun = DolfinFunctionWrapper(self.model.mesh.get_function_space())
#        vec = fun.vector()
#        (nd,) = numpy.shape(data)
#        if nd == len(vec):
#            for i in range(nd):
#                vec[i]=data[i]
#        else:
#            #v2d= self.get_v2d()
#            for i in range(len(vec)):
#                vec[i] = data[i] # shouldn't we use v2d or d2v here?  But it doesn't work if I do.
#        fun.display(opacity=opacity, wireframe=wireframe, width=width, camera=camera)
        raise Exception("TODO")


class ResultError(Exception):
    pass

