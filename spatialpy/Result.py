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


    def __getstate__(self):
        """ Used by pickle to get state when pickling. We need to read the contents of the
        output file since we can't pickle file objects. """
        #TODO
        raise Exception('TODO: spatialpy.Result.__getstate__()');

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """
        #TODO
        raise Exception('TODO: spatialpy.Result.__setstate__()');


    def read_step(self, step_num, debug=False):
        """ Read the data for simulation step 'step_num'. """
        reader = VTKReader(debug=debug)
        filename = os.path.join(self.result_dir, "output{0}.vtk".format(step_num))
        #print("read_step({0}) opening '{1}'".format(step_num, filename))
        reader.setfilename(filename)
        reader.readfile()
        if reader.getpoints() is None or reader.getarrays() is None:
            raise ResultError("read_step(step_num={0}): got data = None".format(step_num))
        points = reader.getpoints()
        vtk_data = reader.getarrays()
        return (points, vtk_data)


    def get_timespan(self):
        self.tspan = numpy.linspace(0,self.model.num_timesteps,
                num=self.model.num_timesteps+1) * self.model.timestep_size
        return self.tspan

    def get_species(self, species, timepoints=None, concentration=False, deterministic=False, debug=False):
        """ Get the populations/concentration values for a given species in the model for 
            one or all timepoints.  
            
            If 'timepoints' is None (default), a matrix of dimension:
            (number of timepoints) x (number of voxels) is returned.  If an integer value is
            given, that value is used to index into the timespan, and that time point is returned
            as a 1D array with size (number of voxel). 

            If concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.

            If deterministic is True, show results for determinstic (instead of stochastic) values
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
            (_, step) = self.read_step(t_ndx, debug=debug)
            if deterministic: 
                ret[ndx,:] = step['C['+spec_name+']']
            elif concentration:
                # concentration = (copy_number/volume)
                # volume = (mass/density)
                ret[ndx,:] = step['D['+spec_name+']'] / (step['mass'] / step['rho'] )
            else:
                ret[ndx,:] = step['D['+spec_name+']']
        if ret.shape[0] == 1:
            ret = ret.flatten()
        return ret

    def get_property(self, property_name, timepoints=None):
        """ Get the property values for a given species in the model for 
            one or all timepoints.  
            
            If 'timepoints' is None (default), a matrix of dimension:
            (number of timepoints) x (number of voxels) is returned.  If an integer value is
            given, that value is used to index into the timespan, and that time point is returned
            as a 1D array with size (number of voxel). 
        """

        t_index_arr = numpy.linspace(0,self.model.num_timesteps, 
                            num=self.model.num_timesteps+1, dtype=int)
        num_voxel = self.model.mesh.get_num_voxels()

        if timepoints is not None:
            if isinstance(timepoints,float):
                raise ResultError("timepoints argument must be an integer, the index of time timespan")
            t_index_arr = t_index_arr[timepoints]
        try:
            num_timepoints = len(t_index_arr)
        except Exception as e:
            t_index_arr = [t_index_arr]
            num_timepoints = 1

        ret = numpy.zeros( (num_timepoints, num_voxel))
        for ndx, t_ndx in enumerate(t_index_arr):
            (_, step) = self.read_step(t_ndx)
            ret[ndx,:] = step[property_name]
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
        #TODO: Check if this still works
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
        #TODO
        raise Exception("todo")




    def display(self, species, time_index, opacity=1.0, wireframe=True, width=500, camera=[0,0,1]):
        """ Plot the trajectory as a PDE style plot. """
        raise Exception("Deprecated")


class ResultError(Exception):
    pass

