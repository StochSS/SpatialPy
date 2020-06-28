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



common_rgb_values=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f',
                   '#bcbd22','#17becf','#ff0000','#00ff00','#0000ff','#ffff00','#00ffff','#ff00ff',
                   '#800000','#808000','#008000','#800080','#008080','#000080','#ff9999','#ffcc99',
                   '#ccff99','#cc99ff','#ffccff','#62666a','#8896bb','#77a096','#9d5a6c','#9d5a6c',
                   '#eabc75','#ff9600','#885300','#9172ad','#a1b9c4','#18749b','#dadecf','#c5b8a8',
                   '#000117','#13a8fe','#cf0060','#04354b','#0297a0','#037665','#eed284','#442244',
                   '#ffddee','#702afb']

common_color_scales = ["Plotly3","Jet","Blues","YlOrRd","PuRd","BuGn","YlOrBr","PuBuGn","BuPu","YlGnBu",
                       "PuBu","GnBu","YlGn","Greens","Reds","Greys","RdPu","OrRd","Purples","Oranges"]


def _plotly_iterate(subdomains, property_name=None, colormap=None):
    import plotly.graph_objs as go

    trace_list = []
    for i, (name, sub_data) in enumerate(subdomains.items()):
        # get point data for trace
        x_data = list(map(lambda point: point[0], sub_data["points"]))
        y_data = list(map(lambda point: point[1], sub_data["points"]))
        z_data = list(map(lambda point: point[2], sub_data["points"]))

        if property_name is not None and property_name == "type":
            marker = {"size":5, "color":common_rgb_values[i]}
        else:
            marker = {"size":5, "color":sub_data["data"], "colorscale":colormap, 
                        "colorbar":{'thickness':20,'title':name}}
        trace = go.Scatter3d(x=x_data, y=y_data, z=z_data, name=name, mode="markers", marker=marker)
        trace_list.append(trace)
    return trace_list


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

    def plot_species(self, species, t_ndx, title=None, concentration=False, deterministic=False, return_plotly_figure=False, width=500, height=500, colormap=None):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

            If concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.

            If deterministic is True, show results for determinstic (instead of stochastic) values

        Attributes
        ----------
        species : str
            A string describing the species to be plotted.
        t_ndx : int
            The time index of the results to be plotted
        title : str
            The title of the graph
        concentration : bool
            Whether or not to plot the data as stochastic concentration, ignored if deterministic is 
            set to True
        deterministic : bool
            Whether or not to plot the data as deterministic
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        colormap : str
            colormap to use
        """
        from plotly.offline import init_notebook_mode, iplot
        
        # read data at time point
        points, data = self.read_step(t_ndx)
        
        # map data to subdomains
        subdomains = {}
        for i, val in enumerate(data['type']):
            name = species
            if deterministic:
                spec_data = data["C[{}]".format(species)][i]
            elif concentration:
                spec_data = data["D[{}]".format(species)][i] / (data['mass'][i] / data['rho'][i])
            else:
                spec_data = data["D[{}]".format(species)][i]

            if name in subdomains.keys():
                subdomains[name]['points'].append(points[i])
                subdomains[name]['data'].append(spec_data)
            else:
                subdomains[name] = {"points":[points[i]], "data":[spec_data]}

        trace_list = _plotly_iterate(subdomains, colormap=colormap)
        
        scene_x = self.model.mesh.xlim[0]/2.5
        scene_y = self.model.mesh.ylim[0]/2.5
        scene_z = self.model.mesh.zlim[0]/2.5
        scene = {"aspectratio": {"x":scene_x,"y":scene_y,"z":scene_y}}
        layout = {"width": width, "height": height, "scene":scene}
        if title is not None:
            layout["title"] = title

        fig = {"data":trace_list, "layout":layout}

        # function for 3D animations
        '''if len(output_data) > 1:
            fig["layout"]["updatemenus"] = [
                {"buttons": [
                    {"args": [None, {"frame": {"duration": 500, "redraw": False},
                                     "fromcurrent": True,
                                     "transition": {"duration": 300, "easing": "quadratic-in-out"}}],
                     "label": "Play",
                     "method": "animate"
                    },
                    {"args": [[None], {"frame": {"duration": 0, "redraw": False},
                                       "mode": "immediate",
                                       "transition": {"duration": 0}}],
                     "label": "Pause",
                     "method": "animate"
                    }],
                 "direction": "left",
                 "pad": {"r": 10, "t": 87},
                 "showactive": False,
                 "type": "buttons",
                 "x": 0.1,
                 "xanchor": "right",
                 "y": 0,
                 "yanchor": "top"
                }]
            
            sliders_dict = {
                "active": 0,
                "yanchor": "top",
                "xanchor": "left",
                "currentvalue": {
                    "font": {"size": 20},
                    "prefix": "Time:",
                    "visible": True,
                    "xanchor": "right"
                },
                "transition": {"duration": 300, "easing": "cubic-in-out"},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": []}
            
            frames = []
            for i in range(0, len(output_data), speed):
                output = output_data[i]
                data = output['fields'][species]
                marker = {"size":5, "color":data, "colorscale":color_scales[species]}
                trace = go.Scatter3d(x=x, y=y, z=z, name=species, mode="markers", marker=marker)
                frame = {"data":[trace], "name":str(i)}
                frames.append(frame)
                
                slider_step = {"args": [[str(i)],
                                        {"frame": {"duration": 100, "redraw": True},
                                         "mode": "immediate",
                                         "transition": {"duration": 100}
                                        }],
                               "label": str(i),
                               "method": "animate"}
                
                sliders_dict['steps'].append(slider_step)
                
            fig["layout"]["sliders"] = [sliders_dict]
            fig["frames"] = frames'''

        if return_plotly_figure:
            return fig
        else:
            iplot(fig)
          
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

    def plot_property(self, property_name, t_ndx, title=None, return_plotly_figure=False, width=500,
                      height=500, colormap=None):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

            If concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.

            If deterministic is True, show results for determinstic (instead of stochastic) values

        Attributes
        ----------
        property_name : str
            A string describing the property to be plotted.
        t_ndx : int
            The time index of the results to be plotted
        title : str
            The title of the graph
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        width : int
            The width of the output window
        height : int
            The height of the output window
        colormap: str
            The style of color maps to use
        """
        from plotly.offline import init_notebook_mode, iplot
        import plotly.graph_objs as go

        # read data at time point
        points, data = self.read_step(t_ndx)

        subdomains = {}
        if property_name == 'type':
            for i, val in enumerate(data['type']):
                name = "type {}".format(val)
                
                if name in subdomains.keys():
                    subdomains[name]['points'].append(points[i])
                    subdomains[name]['data'].append(data[property_name][i])
                else:
                    subdomains[name] = {"points":[points[i]], "data":[data[property_name][i]]}
        else:
            subdomains[property_name] = {
                "points": points,
                "data" : data[property_name]
            }

        trace_list = _plotly_iterate(subdomains, property_name=property_name, colormap=colormap)

        scene_x = self.model.mesh.xlim[0]/2.5
        scene_y = self.model.mesh.ylim[0]/2.5
        scene_z = self.model.mesh.zlim[0]/2.5
        scene = {"aspectratio": {"x":scene_x,"y":scene_y,"z":scene_y}}
        layout = {"width": width, "height": width, "scene":scene}
        if title is not None:
            layout["title"] = title

        fig = {"data":trace_list, "layout":layout}

        if return_plotly_figure:
            return fig
        else:
            iplot(fig)

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

