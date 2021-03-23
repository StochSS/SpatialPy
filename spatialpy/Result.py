import filecmp
import math
import os
import pickle
import shutil
import subprocess

import numpy

from spatialpy.Model import *
from spatialpy.VTKReader import VTKReader

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

    def __init__(self, model=None, result_dir=None, loaddata=False):
        self.model = model
        self.U = None
        self.tspan = None
        self.data_is_loaded = False
        self.success = False
        self.stdout = None
        self.stderr = None
        self.timeout = False
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


    def __eq__(self, other):
        """ Compare Result object's output for equality. This does _NOT_ compare objects themselves

        Params:
            self: Results object
            other: Results object to compare against
        Return:
            bool """

        if isinstance(other, Result) and self.result_dir != None and other.result_dir != None:
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
        """ Compare Result object's output for inequality. This does _NOT_ compare objects themselves.
            This inverts the logic in __eq__().

        Params:
            self: Results object
            other: Results object to compare against
        Return:
            bool """

        return not self.__eq__(other)

    def __getstate__(self):
            """ Used by pickle to get state when pickling. """

            state = {}
            for key, item in self.__dict__.items():
                resultdict = OrderedDict()
                # Pickle all Result output files
                # This does not perserve file metadata like permissions.
                # In the future we should probably at least perserve permissions.
                try:
                    for root, _, file in os.walk(self.result_dir):
                        for filename in file:
                            fd = open(os.path.join(root, filename), 'rb')
                            fd.seek(0)
                            resultdict[filename] = fd.read()
                            fd.close()
                    state['results_output'] = resultdict
                except Exception as e:
                    raise Exception("Error pickling model, could not pickle the Result output files: "+str(e))
                state[key] = item

            return state

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """

        self.__dict__ = state

        # Recreate the Result output files
        # This does not restore file metadata like permissions.
        # In the future we should probably at least restore permissions.
        try:
            results_output = state['results_output']
            if not os.path.exists(state['result_dir']):
                os.mkdir(state['result_dir'])
                for filename, contents in results_output.items():
                    fd = open(os.path.join(state['result_dir'], filename), 'wb')
                    fd.seek(0)
                    fd.write(contents)
                    fd.close()
        except Exception as e:
            raise Exception("Error unpickling model, could not recreate the Result output files: "+str(e))

    def read_step(self, step_num, debug=False):
        """ Read the data for simulation step 'step_num'. """
        reader = VTKReader(debug=debug)
        num = int(step_num * self.model.output_freq)
        filename = os.path.join(self.result_dir, "output{0}.vtk".format(num))
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
                num=math.ceil(self.model.num_timesteps/self.model.output_freq)+1) * self.model.timestep_size
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
        num_voxel = self.model.domain.get_num_voxels()

        if isinstance(species,str):
            spec_name = species
        else:
            spec_name = species.name

        if spec_name not in self.model.listOfSpecies.keys():
            raise ResultError("Species '{0}' not found".format(spec_name))

        #t_index_arr = numpy.linspace(0,self.model.num_timesteps,
        #                    num=self.model.num_timesteps+1, dtype=int)
        t_index_arr = self.get_timespan()

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

    def plot_species(self, species, t_ndx=0, concentration=False, deterministic=False, width=500, height=500, colormap=None, size=5, title=None,
                     animated=False, t_ndx_list=None, speed=1, f_duration=500, t_duration=300, return_plotly_figure=False,
                     use_matplotlib=False, mpl_width=6.4, mpl_height=4.8,
                     debug=False):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

            If concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.

            If deterministic is True, show results for determinstic (instead of stochastic) values

        Attributes
        ----------
        species : str
            A string describing the species to be plotted.
        t_ndx : int
            The time index of the results to be plotted, ignored if animated is set to True
        concentration : bool
            Whether or not to plot the data as stochastic concentration, ignored if deterministic is
            set to True
        deterministic : bool
            Whether or not to plot the data as deterministic
        width: int (default 500)
            Width in pixels of output plot box
        height: int (default 500)
            Height in pixels of output plot box
        colormap : str
            colormap to use.  Plotly specification, valid values: "Plotly3","Jet","Blues","YlOrRd",
                "PuRd","BuGn","YlOrBr","PuBuGn","BuPu","YlGnBu", "PuBu","GnBu","YlGn","Greens","Reds",
                "Greys","RdPu","OrRd","Purples","Oranges".
        size : int
            Size in pixels of the particle
        title : str
            The title of the graph
        animated : bool
            Whether or not the plot is a 3D animation, ignored if use_matplotlib True
        t_ndx_list : list
            The list of time indeces of the results to be plotted, ignored if animated is
            False (default)
        speed : int
            The interval of the time indeces of the results to be plotted (animated plots only)
        f_duration : int
            The duration of time that a frame is displayed
        t_duration : int
            The duration of time to execute the transition between frames
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        use_matplotlib : bool
            whether or not to plot the proprties results using matplotlib.
        mpl_width: int (default 6.4)
            Width in inches of output plot box
        mpl_height: int (default 4.8)
            Height in inches of output plot box
        debug: bool
            output debugging info
        """
        from plotly.offline import init_notebook_mode, iplot

        if(t_ndx < 0):
            t_ndx = len(self.get_timespan()) + t_ndx

        if animated and t_ndx_list is None:
            t_ndx_list = list(range(int(self.model.num_timesteps/self.model.output_freq)+1))

        spec_name = "C[{0}]".format(species) if deterministic else "D[{0}]".format(species)

        # read data at time point
        time_index = t_ndx_list[0] if animated else t_ndx
        points, data = self.read_step(time_index, debug=debug)

        if use_matplotlib:
            import matplotlib.pyplot as plt

            if (deterministic or not concentration):
                d = data[spec_name]
            else:
                d = data[spec_name] / (data['mass'] / data['rho'])
            if colormap is None:
                colormap = "viridis"

            plt.figure(figsize=(mpl_width,mpl_height))
            plt.scatter(points[:,0],points[:,1],c=d,cmap=colormap)
            plt.axis('scaled')
            plt.colorbar()
            if title is not None:
                plt.title(title)
            plt.grid(linestyle='--', linewidth=1)
            plt.plot()
            return

        # map data to types
        types = {}
        for i, val in enumerate(data['type']):
            name = species
            if deterministic or not concentration:
                spec_data = data[spec_name][i]
            else:
                spec_data = data[spec_name][i] / (data['mass'][i] / data['rho'][i])

            if name in types.keys():
                types[name]['points'].append(points[i])
                types[name]['data'].append(spec_data)
            else:
                types[name] = {"points":[points[i]], "data":[spec_data]}

        is_2d = self.model.domain.zlim[0] == self.model.domain.zlim[1]

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

            _data = data[spec_name] if deterministic or not concentration else data[spec_name] / (data['mass'] / data['rho'])
            cmin = min(_data)
            cmax = max(_data)
            for i in range(1, len(t_ndx_list), speed):
                _, _data = self.read_step(t_ndx_list[i])
                _data = _data[spec_name] if deterministic or not concentration else _data[spec_name] / (_data['mass'] / _data['rho'])
                if min(_data) - 0.1 < cmin:
                    cmin = min(_data) - 0.1
                if max(_data) + 0.1 > cmax:
                    cmax = max(_data) + 0.1

            frames = []
            for index in range(0, len(t_ndx_list), speed):
                points, data = self.read_step(t_ndx_list[index])

                # map data to types
                types = {}
                for i, val in enumerate(data['type']):
                    name = species
                    if deterministic or not concentration:
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
        """ Get the property values for a given species in the model for
            one or all timepoints.

            If 'timepoints' is None (default), a matrix of dimension:
            (number of timepoints) x (number of voxels) is returned.  If an integer value is
            given, that value is used to index into the timespan, and that time point is returned
            as a 1D array with size (number of voxel).
        """

        t_index_arr = numpy.linspace(0,self.model.num_timesteps,
                            num=self.model.num_timesteps+1, dtype=int)
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
            (_, step) = self.read_step(t_ndx)
            ret[ndx,:] = step[property_name]
        if ret.shape[0] == 1:
            ret = ret.flatten()
        return ret

    def plot_property(self, property_name, t_ndx=0, p_ndx=0, width=500, height=500, colormap=None, size=5, title=None,
                      animated=False, t_ndx_list=None, speed=1, f_duration=500, t_duration=300, return_plotly_figure=False,
                      use_matplotlib=False, mpl_width=6.4, mpl_height=4.8):
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
        p_ndx : int
            The property index of the results to be plotted
        width: int (default 500)
            Width in pixels of output plot box
        height: int (default 500)
            Height in pixels of output plot box
        colormap : str
            colormap to use.  Plotly specification, valid values: "Plotly3","Jet","Blues","YlOrRd",
                "PuRd","BuGn","YlOrBr","PuBuGn","BuPu","YlGnBu", "PuBu","GnBu","YlGn","Greens","Reds",
                "Greys","RdPu","OrRd","Purples","Oranges".
        size : int
            Size in pixels of the particle
        title : str
            The title of the graph
        animated : bool
            Whether or not the plot is a 3D animation, ignored if use_matplotlib True
        t_ndx_list : list
            The list of time indeces of the results to be plotted, ignored if animated is
            False (default)
        speed : int
            The interval of the time indeces of the results to be plotted (animated plots only)
        f_duration : int
            The duration of time that a frame is displayed
        t_duration : int
            The duration of time to execute the transition between frames
        return_plotly_figure : bool
            whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        use_matplotlib : bool
            whether or not to plot the proprties results using matplotlib.
        mpl_width: int (default 6.4)
            Width in inches of output plot box
        mpl_height: int (default 4.8)
            Height in inches of output plot box
        """
        if(t_ndx < 0):
            t_ndx = len(self.get_timespan()) + t_ndx

        if animated and t_ndx_list is None:
            t_ndx_list = list(range(int(self.model.num_timesteps/self.model.output_freq)+1))

        # read data at time point
        time_index = t_ndx_list[0] if animated else t_ndx
        points, data = self.read_step(time_index)

        if use_matplotlib:
            import matplotlib.pyplot as plt

            if (property_name == 'v'):
                d = data[property_name]
                d = [d[i][p_ndx] for i in range(0,len(d))]
            else:
                d = data[property_name]
            if colormap is None:
                colormap = "viridis"

            plt.figure(figsize=(mpl_width,mpl_height))
            plt.scatter(points[:,0],points[:,1],c=d,cmap=colormap)
            plt.axis('scaled')
            plt.colorbar()
            if title is not None:
                plt.title(title)
            plt.grid(linestyle='--', linewidth=1)
            plt.plot()
            return

        from plotly.offline import init_notebook_mode, iplot

        types = {}
        if property_name == 'type':
            for i, val in enumerate(data['type']):
                name = "type {}".format(val)

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

        is_2d = self.model.domain.zlim[0] == self.model.domain.zlim[1]

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
                _, _data = self.read_step(t_ndx_list[i])
                _cmin = min(_data[property_name]) if property_name != "v" else min(_data[property_name], key=lambda val: val[p_ndx])[p_ndx]
                if _cmin - 0.1 < cmin:
                    cmin = _cmin - 0.1
                _cmax = max(_data[property_name]) if property_name != "v" else max(_data[property_name], key=lambda val: val[p_ndx])[p_ndx]
                if _cmax + 0.1 > cmax:
                    cmax = _cmax + 0.1

            frames = []
            for index in range(0, len(t_ndx_list), speed):
                points, data = self.read_step(t_ndx_list[index])

                # map data to types
                types = {}
                if property_name == 'type':
                    for i, val in enumerate(data['type']):
                        name = "type {}".format(val)

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
            The columns of mesh.csv are: 'Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Type'.
            The columns of species_S.csv are: 'Time', 'Voxel 0', Voxel 1', ... 'Voxel N'.
        """
        #TODO: Check if this still works
        import csv
        subprocess.call(["mkdir", "-p", folder_name])
        #['Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Type']
        with open(os.path.join(folder_name,'mesh.csv'), 'w+') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['Voxel ID', 'X', 'Y', 'Z', 'Volume', 'Type'])
            vol = self.model.get_solver_datastructure()['vol']
            for ndx in range(self.model.domain.get_num_voxels()):
                row = [ndx]+self.model.domain.coordinates()[ndx,:].tolist()+[vol[ndx]]+[self.model.domain.type[ndx]]
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
