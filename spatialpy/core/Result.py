'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import csv
import filecmp
import math
import os
import shutil
import tempfile
import sys

import numpy

from spatialpy.core.Model import *
from spatialpy.core.VTKReader import VTKReader
from spatialpy.core.spatialpyError import ResultError

try:
    import vtk
except ImportError as e:
    pass

common_rgb_values=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f',
                   '#bcbd22','#17becf','#ff0000','#00ff00','#0000ff','#ffff00','#00ffff','#ff00ff',
                   '#800000','#808000','#008000','#800080','#008080','#000080','#ff9999','#ffcc99',
                   '#ccff99','#cc99ff','#ffccff','#62666a','#8896bb','#77a096','#9d5a6c','#9d5a6c',
                   '#eabc75','#ff9600','#885300','#9172ad','#a1b9c4','#18749b','#dadecf','#c5b8a8',
                   '#000117','#13a8fe','#cf0060','#04354b','#0297a0','#037665','#eed284','#442244',
                   '#ffddee','#702afb']

common_color_scales = ["Plotly3","Jet","Blues","YlOrRd","PuRd","BuGn","YlOrBr","PuBuGn","BuPu","YlGnBu",
                       "PuBu","GnBu","YlGn","Greens","Reds","Greys","RdPu","OrRd","Purples","Oranges"]


def _plotly_iterate(types, size=5, property_name=None, cmin=None, cmax=None, colormap=None, is_2d=False):
    """ Helper method used by plotly display methods. """
    import plotly.graph_objs as go

    trace_list = []
    for i, (name, sub_data) in enumerate(types.items()):
        # get point data for trace
        x_data = list(map(lambda point: point[0], sub_data["points"]))
        y_data = list(map(lambda point: point[1], sub_data["points"]))
        z_data = list(map(lambda point: point[2], sub_data["points"]))

        if property_name is not None and property_name == "type":
            marker = {"size":size, "color":common_rgb_values[i]}
        else:
            if colormap is None:
                colormap = common_color_scales[i]
            marker = {"size":size, "color":sub_data["data"], "colorscale":colormap,
                        "colorbar":{'thickness':20,'title':name}}
            if cmin is not None and cmax is not None:
                marker["cmin"] = cmin
                marker["cmax"] = cmax

        if is_2d:
            trace = go.Scatter(x=x_data, y=y_data, name=name, mode="markers", marker=marker)
        else:
            trace = go.Scatter3d(x=x_data, y=y_data, z=z_data, name=name, mode="markers", marker=marker)
        trace_list.append(trace)
    return trace_list

class Result():
    """ Result object for a URDME simulation. """

    def __init__(self, model=None, result_dir=None):
        self.model = model
        self.tspan = None
        self.success = False
        self.stdout = None
        self.stderr = None
        self.timeout = False
        self.official_vtk = False
        self.result_dir = result_dir


    def __eq__(self, other):
        """ Compare Result object's simulation data for equality. This does _NOT_ compare objects themselves.

        Attributes
        ----------
            other: Result
                Results object to compare against

        Return
        ----------
            bool:
                Whether or not the Results object's simulation data is equal
        """

        if isinstance(other, Result) and self.result_dir and other.result_dir:
            # Compare contents, not shallow compare
            filecmp.cmpfiles.__defaults__ = (False,)
            dircmp = filecmp.dircmp(self.result_dir, other.result_dir)
            # Raise exception if funny_files
            assert not dircmp.funny_files
            if not (dircmp.left_only or dircmp.right_only or dircmp.funny_files or dircmp.diff_files):
                return True
            return False
        return NotImplemented

    def __ne__(self, other):
        """ Compare Result object's simulation data for inequality. This does _NOT_ compare objects themselves.

        Attributes
        ----------
            other: Result
                Results object to compare against

        Return
        ----------
            bool:
                Whether or not the Results object's simulation data is unequal
        """

        return not self.__eq__(other)

    def __getstate__(self):
        """ Used by pickle to get state when pickling. """

        state = {}
        for key, item in self.__dict__.items():
            resultdict = OrderedDict()

            try:
                for root, _, file in os.walk(self.result_dir):
                    for filename in file:
                        with open(os.path.join(root, filename), 'rb') as fd:
                            fd.seek(0)
                            resultdict[filename] = fd.read()
                state['results_output'] = resultdict
            except Exception as e:
                raise ResultError("Error pickling model, could not pickle the Result output files: "+str(e))
            state[key] = item

        return state

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """

        self.__dict__ = state

        try:
            results_output = state['results_output']
            state['result_dir'] = tempfile.mkdtemp(
                prefix='spatialpy_result_', dir=os.environ.get('SPATIALPY_TMPDIR'))

            for filename, contents in results_output.items():
                with open(os.path.join(state['result_dir'], filename), 'wb') as fd:
                    fd.seek(0)
                    fd.write(contents)
        except Exception as e:
            raise ResultError("Error unpickling model, could not recreate the Result output files: "+str(e))

    def __del__(self):
        """ Deconstructor. """
        try:
            if self.result_dir is not None:
                try:
                    shutil.rmtree(self.result_dir)
                except OSError as e:
                    print("Could not delete '{0}'".format(self.result_dir))
        except Exception as e:
            pass

    def read_step(self, step_num, debug=False):
        """ Read the data for simulation step 'step_num'. Returns a tuple containing a numpy.ndarray \
            of point coordinates [0] along with a dictionary of property and species data [1]

        :param step_num: The index in the timespan 
        :type step_num: int
        :param debug: Whether or not debug information should be printed
        :type debug: bool (default False)

        :rtype: tuble

        """

        if debug:
            print("read_step({0}) ".format(step_num), end='')
        filename = os.path.join(self.result_dir, "output{0}.vtk".format(step_num))

        if debug:
            print("opening '{0}'".format(filename))

        if self.official_vtk:
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(filename)
            reader.Update()
            data = reader.GetOutput()

            if data is not None:
                points = numpy.array(data.GetPoints().GetData())
                pd = data.GetPointData()
                vtk_data = {}

                for i in range(pd.GetNumberOfArrays()):
                    if pd.GetArrayName(i) is None:
                        break

                    if debug:
                        print(i,pd.GetArrayName(i))

                    vtk_data[pd.GetArrayName(i)] = numpy.array(pd.GetArray(i))
        else:
            reader = VTKReader(filename=filename, debug=debug)
            reader.readfile()
            points = reader.getpoints()
            vtk_data = reader.getarrays()

        if points is None or vtk_data is None:
            raise ResultError("read_step(step_num={0}): got data = None".format(step_num))

        return (points, vtk_data)

    def get_timespan(self):
        """  Get the model time span. Returns a numpy array containing the time span of the model

            :rtype: numpy.ndarray:

        """

        self.tspan = self.model.tspan
        return self.tspan

    def get_species(self, species, timepoints=None, concentration=False, deterministic=False, debug=False):
        """ Get the populations/concentration values for a given species in the model for
            one or all timepoints. Returns a numpy array containing the species population/concentration values

                :param species: A species in string or dictionary form to retreive information about
                :type species: str | dict
                :param timepoints: A time point where the information should be retreived from.
                    If 'timepoints' is None (default), a matrix of dimension:
                    (number of timepoints) x (number of voxels) is returned.
                    If an integer value is given, that value is used to index into the timespan, and that time point is returned
                    as a 1D array with size (number of voxel).
                :type timepoints: int (default None)
                :param concentration: Whether or not the species is a concentration (True) or population (False)
                    If concentration is False (default), the integer, raw, trajectory data is returned.
                    If set to True, the concentration (=copy_number/volume) is returned.
                :type concentration: bool (default False)
                :param deterministic: Whether or not the species is deterministic (True) or stochastic (False)
                :type deterministic: bool (default False)
                :param debug: Whether or not debug information should be printed
                :type debug: bool (default False)

                :rtype: numpy.ndarray
        """

        species_map = self.model.species_map
        num_species = self.model.get_num_species()
        num_voxel = self.model.domain.get_num_voxels()

        if isinstance(species,str):
            spec_name = species
        else:
            spec_name = species.name

        if spec_name not in self.model.listOfSpecies.keys():
            raise ResultError("Species '{0}' not found".format(spec_name))

        lt=len(self.get_timespan())-1
        t_index_arr = numpy.linspace(0,lt,num=lt+1,dtype=int)
        
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

    def plot_species(self, species, t_ndx=None, t_val=None, concentration=False, 
                     deterministic=False, width=None, height=None, colormap=None,
                     size=5, title=None, animated=False, t_ndx_list=None, speed=1,
                     f_duration=500, t_duration=300, return_plotly_figure=False,
                     use_matplotlib=False, debug=False):
        """ Plots the Results using plotly or matplotlib. Can be viewed in a Jupyter Notebook.

            If concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.

            If deterministic is True, show results for determinstic (instead of stochastic) values

        :param species: A string describing the species to be plotted.
        :type species: str
        :param t_ndx: The time index of the results to be plotted, ignored if animated is set to True
        :type t_ndx: int
        :param t_val: The time value of the results to be plotted, ignored if animated is set to True
        :type t_val: float
        :param concentration: Whether or not to plot the data as stochastic concentration, ignored if deterministic is
            set to True
        :type concentration: bool
        :param deterministic: Whether or not to plot the data as deterministic
        :type deterministic: bool
        :param width: Width in pixels of output plot box
        :type width: int (default 500)
        :param height: Height in pixels of output plot box
        :type height: int (default 500)
        :param colormap: colormap to use.  Plotly specification, valid values: "Plotly3","Jet","Blues","YlOrRd",
                "PuRd","BuGn","YlOrBr","PuBuGn","BuPu","YlGnBu", "PuBu","GnBu","YlGn","Greens","Reds",
                "Greys","RdPu","OrRd","Purples","Oranges".
        :type colormap: str
        :param size: Size in pixels of the particle
        :type size: int
        :param title: The title of the graph
        :type title: str
        :param animated: Whether or not the plot is a 3D animation, ignored if use_matplotlib True
        :type animated: bool
        :param t_ndx_list: The list of time indeces of the results to be plotted, ignored if animated is
            False (default)
        :type t_ndx_list: list
        :param speed: The interval of the time indeces of the results to be plotted (animated plots only)
        :type speed: int
        :param f_duration: The duration of time that a frame is displayed
        :type f_duration: int
        :param t_duration: The duration of time to execute the transition between frames
        :type t_duration: int
        :param return_plotly_figure: whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        :type return_plotly_figure: bool
        :param use_matplotlib: whether or not to plot the proprties results using matplotlib.
        :type use_matplotlib: bool
        :param debug: output debugging info
        :type debug: bool
        """

        from plotly.offline import init_notebook_mode, iplot


        time_index_list = self.get_timespan()

        if animated and t_ndx_list is None:
            t_ndx_list = time_index_list 

        spec_name = "C[{0}]".format(species) if deterministic else "D[{0}]".format(species)

        # read data at time point
        if animated:
            t_ndx = 0
        elif t_val is None and t_ndx is None:
            t_ndx = 0 # default value
        elif t_val is not None and t_ndx is not None:
            raise ResultError("t_ndx and t_val can not both be set.")
        elif t_val is not None:
            if t_val in time_index_list:
                for n,t in enumerate(time_index_list):
                    if t == t_val:
                        t_ndx = n
                        break
            else:
                raise ResultError("time value (t_val) value given is not a valid result timepoint") 
        else:
            if t_ndx != int(t_ndx):
                raise ResultError("t_ndx must be an integer")
            if t_ndx < 0:
                t_ndx = len(self.get_timespan()) + t_ndx
        points, data = self.read_step(t_ndx, debug=debug)

        if use_matplotlib:
            import matplotlib.pyplot as plt

            if width is None:
                width = 6.4
            if height is None:
                height = 4.8

            if deterministic or not concentration:
                d = data[spec_name]
            else:
                d = data[spec_name] / (data['mass'] / data['rho'])
            if colormap is None:
                colormap = "viridis"

            plt.figure(figsize=(width, height))
            plt.scatter(points[:,0],points[:,1],c=d,cmap=colormap)
            plt.axis('scaled')
            plt.colorbar()
            if title is not None:
                plt.title(title)
            plt.grid(linestyle='--', linewidth=1)
            plt.plot()
            return

        if width is None:
            width = 500
        if height is None:
            height = 500

        # map data to types
        types = {}
        for i, val in enumerate(data['type']):
            name = species
            if (deterministic or not concentration):
                spec_data = data[spec_name][i]
            else:
                spec_data = data[spec_name][i] / (data['mass'][i] / data['rho'][i])

            if name in types.keys():
                types[name]['points'].append(points[i])
                types[name]['data'].append(spec_data)
            else:
                types[name] = {"points":[points[i]], "data":[spec_data]}

        is_2d = self.model.domain.dimensions == 2

        trace_list = _plotly_iterate(types, size=size, colormap=colormap, is_2d=is_2d)

        scene = {
            "aspectmode": 'data',
        }
        layout = {"width": width, "height": width, "scene":scene,
                  "xaxis":{"range":self.model.domain.xlim}, "yaxis":{"range":self.model.domain.ylim}
                 }
        if title is not None:
            layout["title"] = title

        fig = {"data":trace_list, "layout":layout}

        # function for 3D animations
        if animated and len(t_ndx_list) > 1:
            fig["layout"]["updatemenus"] = [
                {"buttons": [
                    {"args": [None, {"frame": {"duration": f_duration, "redraw": False},
                                     "fromcurrent": True,
                                     "transition": {"duration": t_duration, "easing": "quadratic-in-out"}}],
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
                "transition": {"duration": t_duration, "easing": "cubic-in-out"},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": []}

            _data = data[spec_name] if (deterministic or not concentration) else data[spec_name] / (data['mass'] / data['rho'])
            cmin = min(_data)
            cmax = max(_data)
            for i in range(1, len(t_ndx_list), speed):
                _, _data = self.read_step(i)
                _data = _data[spec_name] if (deterministic or not concentration) else _data[spec_name] / (_data['mass'] / _data['rho'])
                if min(_data) - 0.1 < cmin:
                    cmin = min(_data) - 0.1
                if max(_data) + 0.1 > cmax:
                    cmax = max(_data) + 0.1

            frames = []
            for index in range(0, len(t_ndx_list), speed):
                points, data = self.read_step(index)

                # map data to types
                types = {}
                for i, val in enumerate(data['type']):
                    name = species
                    if (deterministic or not concentration):
                        spec_data = data[spec_name][i]
                    else:
                        spec_data = data[spec_name][i] / (data['mass'][i] / data['rho'][i])

                    if name in types.keys():
                        types[name]['points'].append(points[i])
                        types[name]['data'].append(spec_data)
                    else:
                        types[name] = {"points":[points[i]], "data":[spec_data]}

                trace_list = _plotly_iterate(types, size=size, colormap=colormap, cmin=cmin, cmax=cmax, is_2d=is_2d)

                frame = {"data":trace_list, "name":str(t_ndx_list[index])}
                frames.append(frame)

                slider_step = {"args": [[str(t_ndx_list[index])],
                                        {"frame": {"duration": f_duration, "redraw": True},
                                         "mode": "immediate",
                                         "transition": {"duration": t_duration}
                                        }],
                               "label": str(t_ndx_list[index]),
                               "method": "animate"}

                sliders_dict['steps'].append(slider_step)

            fig["layout"]["sliders"] = [sliders_dict]
            fig["frames"] = frames

        if return_plotly_figure:
            return fig
        else:
            init_notebook_mode(connected=True)
            iplot(fig)

    def get_property(self, property_name, timepoints=None):
        """ Get the property values for a given species in the model for \
            one or all timepoints. If 'timepoints' is None (default), a matrix of dimension: \
            (number of timepoints) x (number of voxels) is returned.  If an integer value is \
            given, that value is used to index into the timespan, and that time point is returned \
            as a 1D array with size (number of voxel).

            :param property_name: A string describing the property to be returned.
            :type property_name: str
            :param timepoints: timespan index to be returned.  Default is None
            :type timepoints: int

            :rtype: numpy.ndarray
        """

        lt=len(self.get_timespan())-1
        t_index_arr = numpy.linspace(0,lt,num=lt+1,dtype=int)
        
        num_voxel = self.model.domain.get_num_voxels()

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
            (_, step) = self.read_step(ndx)
            ret[ndx,:] = step[property_name]
        if ret.shape[0] == 1:
            ret = ret.flatten()
        return ret

    def plot_property(self, property_name, t_ndx=None, t_val=None, p_ndx=0, width=None, height=None,
                      colormap=None, size=5, title=None, animated=False, t_ndx_list=None, speed=1,
                      f_duration=500, t_duration=300, included_types_list=None, return_plotly_figure=False,
                      use_matplotlib=False, debug=False):
        """
        Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

        :param property_name: A string describing the property to be plotted.
        :type property_name: str
        :param t_ndx: The time index of the results to be plotted
        :type t_ndx: int
        :param t_val: The time value of the results to be plotted, ignored if animated is set to True
        :type t_val: float
        :param p_ndx: The property index of the results to be plotted
        :type p_ndx: int
        :param width: Width in pixels of output plot box or for matplotlib inches of output plot box
        :type width: int (default 500)
        :param height: Height in pixels of output plot box or for matplotlib inches of output plot box
        :type height: int (default 500)
        :param colormap: colormap to use.  Plotly specification, valid values: "Plotly3","Jet","Blues","YlOrRd",
                "PuRd","BuGn","YlOrBr","PuBuGn","BuPu","YlGnBu", "PuBu","GnBu","YlGn","Greens","Reds",
                "Greys","RdPu","OrRd","Purples","Oranges".
        :type colormap: str
        :param size: Size in pixels of the particle
        :type size: int
        :param title: The title of the graph
        :type title: str
        :param animated: Whether or not the plot is a 3D animation, ignored if use_matplotlib True
        :type animated: bool
        :param t_ndx_list: The list of time indeces of the results to be plotted, ignored if animated is
            False (default)
        :type t_ndx_list: list
        :param speed: The interval of the time indeces of the results to be plotted (animated plots only)
        :type speed: int
        :param f_duration: The duration of time that a frame is displayed
        :type f_duration: int
        :param t_duration: The duration of time to execute the transition between frames
        :type t_duration: int
        :param included_types_list: A list of ints describing which types to include. By default displays all types.
        :type included_types_list: list
        :param return_plotly_figure: whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        :type return_plotly_figure: bool
        :param use_matplotlib: whether or not to plot the proprties results using matplotlib.
        :type use_matplotlib: bool
        :param debug: output debugging info
        :type debug: bool
        """

        time_index_list = self.get_timespan()

        if animated and t_ndx_list is None:
            t_ndx_list = time_index_list

        # read data at time point
        if animated:
            t_ndx = 0
        elif t_val is None and t_ndx is None:
            t_ndx = 0 # default value
        elif t_val is not None and t_ndx is not None:
            raise ResultError("t_ndx and t_val can not both be set.")
        elif t_val is not None:
            if t_val in time_index_list:
                for n,t in enumerate(time_index_list):
                    if t == t_val:
                        t_ndx = n
                        break
            else:
                raise ResultError("time value (t_val) value given is not a valid result timepoint")
        else:
            if t_ndx != int(t_ndx):
                raise ResultError("t_ndx must be an integer")
            if t_ndx < 0:
                t_ndx = len(self.get_timespan()) + t_ndx
        points, data = self.read_step(t_ndx, debug=debug)

        if use_matplotlib:
            import matplotlib.pyplot as plt

            if width is None:
                width = 6.4
            if height is None:
                height = 4.8

            if property_name == "type" and included_types_list is not None:
                coords = []
                d = []
                for i, val in enumerate(data[property_name]):
                    if val in included_types_list:
                        coords.append(points[i])
                        d.append(val)
                points = numpy.array(coords)
            elif property_name == 'v':
                d = data[property_name]
                d = [d[i][p_ndx] for i in range(0,len(d))]
            else:
                d = data[property_name]
            if colormap is None:
                colormap = "viridis"

            plt.figure(figsize=(width, height))
            plt.scatter(points[:,0],points[:,1],c=d,cmap=colormap)
            plt.axis('scaled')
            plt.colorbar()
            if title is not None:
                plt.title(title)
            plt.grid(linestyle='--', linewidth=1)
            plt.plot()
            return

        if width is None:
            width = 500
        if height is None:
            height = 500

        from plotly.offline import init_notebook_mode, iplot

        types = {}
        if property_name == 'type':
            for i, val in enumerate(data['type']):
                name = "type {}".format(val)

                if included_types_list is None or val in included_types_list:
                    if name in types.keys():
                        types[name]['points'].append(points[i])
                        types[name]['data'].append(data[property_name][i])
                    else:
                        types[name] = {"points":[points[i]], "data":[data[property_name][i]]}
        elif property_name == 'v':
            types[property_name] = {
                "points": points,
                "data" : [data[property_name][i][p_ndx] for i in range(0,len(data[property_name]))]
            }
        else:
            types[property_name] = {
                "points": points,
                "data" : data[property_name]
            }

        is_2d = self.model.domain.dimensions == 2

        trace_list = _plotly_iterate(types, size=size, property_name=property_name,
                                     colormap=colormap, is_2d=is_2d)

        scene = {
            "aspectmode": 'data',
        }
        layout = {"width": width, "height": width, "scene":scene,
                  "xaxis":{"range":self.model.domain.xlim}, "yaxis":{"range":self.model.domain.ylim}
                 }

        if title is not None:
            layout["title"] = title

        fig = {"data":trace_list, "layout":layout}

        # function for 3D animations
        if animated and len(t_ndx_list) > 1:
            fig["layout"]["updatemenus"] = [
                {"buttons": [
                    {"args": [None, {"frame": {"duration": f_duration, "redraw": False},
                                     "fromcurrent": True,
                                     "transition": {"duration": t_duration, "easing": "quadratic-in-out"}}],
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
                "transition": {"duration": t_duration, "easing": "cubic-in-out"},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": []}

            cmin = min(data[property_name]) if property_name != "v" else min(data[property_name], key=lambda val: val[p_ndx])[p_ndx]
            cmax = max(data[property_name]) if property_name != "v" else max(data[property_name], key=lambda val: val[p_ndx])[p_ndx]
            for i in range(1, len(t_ndx_list), speed):
                _, _data = self.read_step(i)
                _cmin = min(_data[property_name]) if property_name != "v" else min(_data[property_name], key=lambda val: val[p_ndx])[p_ndx]
                if _cmin - 0.1 < cmin:
                    cmin = _cmin - 0.1
                _cmax = max(_data[property_name]) if property_name != "v" else max(_data[property_name], key=lambda val: val[p_ndx])[p_ndx]
                if _cmax + 0.1 > cmax:
                    cmax = _cmax + 0.1

            frames = []
            for index in range(0, len(t_ndx_list), speed):
                points, data = self.read_step(index)

                # map data to types
                types = {}
                if property_name == 'type':
                    for i, val in enumerate(data['type']):
                        name = "type {}".format(val)

                        if included_types_list is None or val in included_types_list:
                            if name in types.keys():
                                types[name]['points'].append(points[i])
                                types[name]['data'].append(data[property_name][i])
                            else:
                                types[name] = {"points":[points[i]], "data":[data[property_name][i]]}
                elif property_name == 'v':
                    types[property_name] = {
                        "points": points,
                        "data" : [data[property_name][i][p_ndx] for i in range(0,len(data[property_name]))]
                    }
                else:
                    types[property_name] = {
                        "points": points,
                        "data" : data[property_name]
                    }

                trace_list = _plotly_iterate(types, size=size, property_name=property_name,
                                     colormap=colormap, cmin=cmin, cmax=cmax, is_2d=is_2d)

                frame = {"data":trace_list, "name":str(t_ndx_list[index])}
                frames.append(frame)

                slider_step = {"args": [[str(t_ndx_list[index])],
                                        {"frame": {"duration": f_duration, "redraw": False},
                                         "mode": "immediate",
                                         "transition": {"duration": t_duration}
                                        }],
                               "label": str(t_ndx_list[index]),
                               "method": "animate"}

                sliders_dict['steps'].append(slider_step)

            fig["layout"]["sliders"] = [sliders_dict]
            fig["frames"] = frames

        if return_plotly_figure:
            return fig
        else:
            init_notebook_mode(connected=True)
            iplot(fig)

    def export_to_csv(self, folder_name=None):
        """ Write the trajectory to a set of CSV files. The first, modelname_mesh.csv, specifies the mesh.
            The other files, modelname_species_S.csv, for species named S, specify trajectory data for each species.
            The columns of modelname_mesh.csv are: 'Voxel ID', 'X', 'Y', 'Z', 'Type', 'Volume', 'Mass', 'Viscosity'
            The columns of modelname_species_S.csv: 'Time', 'Voxel 0', Voxel 1', ... 'Voxel N'.

            :type folder_name: str (default current working directory)
            :param folder_name: A path where the vtk files will be written, created if non-existant.
                If no path is provided current working directory is used.
        """

        if not folder_name:
            folder_name = os.path.abspath(os.getcwd())
        elif not os.path.exists(folder_name):
            os.mkdir(folder_name)

        #['Voxel ID', 'X', 'Y', 'Z', 'Type', 'Volume', 'Mass', 'Viscosity']
        with open(os.path.join(folder_name, self.model.name + '_mesh.csv'), 'w+') as csvfile:
            mesh = self.model.mesh
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['Voxel ID', 'X', 'Y', 'Z', 'Type', 'Volume', 'Mass', 'Viscosity'])
            for ndx in range(len(mesh.vertices)):
                writer.writerow([ndx] + mesh.coordinates()[ndx,:].tolist() + [mesh.type[ndx]] \
                    + [mesh.vol[ndx]] + [mesh.mass[ndx]] + [mesh.nu[ndx]])

        for species in self.model.listOfSpecies:
            #['Voxel', 'Time 0', Time 1', ... 'Time N']
            with open(os.path.join(folder_name, self.model.name + '_species_{0}.csv'.format(species)), 'w+') as csvfile:
                data = self.get_species(species)
                (num_time, num_vox) = data.shape
                writer = csv.writer(csvfile, delimiter=',')
                header_row = ['Voxel']
                for time in range(num_time):
                    header_row.append('Time {0}'.format(time))
                writer.writerow(header_row)
                for voxel in range(num_vox):
                    writer.writerow([voxel] + data[:,voxel].tolist())

    def export_to_vtk(self, timespan, folder_name=None):
        """ Write the trajectory to a collection of vtk files.
            The exported data is #molecules/volume, where the volume unit is implicit from the mesh dimension.
            Not currently implemented."""
            # TODO
        raise ResultError("Not implemented.")