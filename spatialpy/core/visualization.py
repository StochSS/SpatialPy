# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2022 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import math
import numpy

try:
    import matplotlib.pyplot as plt
    mpl_import_err = None
except ImportError as err:
    mpl_import_err = err

from spatialpy.core.spatialpyerror import VisualizationError

common_rgb_values=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f',
                   '#bcbd22','#17becf','#ff0000','#00ff00','#0000ff','#ffff00','#00ffff','#ff00ff',
                   '#800000','#808000','#008000','#800080','#008080','#000080','#ff9999','#ffcc99',
                   '#ccff99','#cc99ff','#ffccff','#62666a','#8896bb','#77a096','#9d5a6c','#9d5a6c',
                   '#eabc75','#ff9600','#885300','#9172ad','#a1b9c4','#18749b','#dadecf','#c5b8a8',
                   '#000117','#13a8fe','#cf0060','#04354b','#0297a0','#037665','#eed284','#442244',
                   '#ffddee','#702afb']

def _get_coords(points):
    labels = ["X-Axis", "Y-Axis", "Z-Axis"]
    coords = []
    axes = []
    i = 0
    while len(coords) < 2 and i < 3:
        vertex = list(map(lambda point: point[i], points))
        if numpy.count_nonzero(numpy.diff(vertex)) > 0:
            coords.append(vertex)
            axes.append(labels[i])
        i += 1
    if len(coords) < 2:
        if labels[0] not in axes:
            coords.insert(0, list(map(lambda point: point[0], points)))
            axes.insert(0, labels[0])
        elif labels[1] not in axes:
            coords.append(list(map(lambda point: point[1], points)))
            axes.append(labels[1])
    return coords, axes

def _validate_mplplot_args(args):
    if args is None:
        return {"figsize": (Visualization.MPL_WIDTH, Visualization.MPL_HEIGHT)}

    supported_splts_args = ["sharex", "sharey", "squeeze"]
    supported_fig_args = ["figsize", "dpi", "facecolor", "edgecolor", "frameon"]

    kwargs = {}
    for name, val in args.items():
        if name in supported_splts_args or name in supported_fig_args:
            kwargs[name] = val
        else:
            from spatialpy.core import log # pylint: disable=import-outside-toplevel
            logmsg = f"Un-supported key word argument: {name} is not currently supported"
            log.warning(logmsg)
    return kwargs

class Visualization():

    MPL_WIDTH = 6.4
    MPL_HEIGHT = 4.8
    MPL_SIZE = 40

    def __init__(self, data):
        self.data = data

    def __get_grid_shape(self, multiple_graphs):
        if isinstance(multiple_graphs, tuple):
            num_subplots = multiple_graphs[0] * multiple_graphs[1]
            if num_subplots < len(self.data.keys()):
                errmsg = f"The shape {multiple_graphs} of the graphs is to small for the given data"
                raise VisualizationError(errmsg)
            nrows = multiple_graphs[0]
            ncols = multiple_graphs[1]
        else:
            nrows = math.ceil(len(self.data.keys()) / 2)
            ncols = 2
        return nrows, ncols

    def __validate_mplscatter_args(self, args, name):
        if args is None:
            return {"s": Visualization.MPL_SIZE}

        supported_args = ["s", "cmap", "color", "marker", "vmin", "vmax"]

        group_args = args[name]
        if "cmap" in group_args and "color" in group_args:
            errmsg = f"scatter args for {name}: 'cmap' and 'color' cannot both be set."
            raise VisualizationError(errmsg)

        kwargs = {}
        for arg_name, val in group_args.items():
            if arg_name in supported_args:
                kwargs[arg_name] = val
            else:
                from spatialpy.core import log # pylint: disable=import-outside-toplevel
                logmsg = f"Un-supported key word argument: {arg_name} is not currently supported"
                log.warning(logmsg)
        if "s" not in kwargs:
            if "size_scale" in self.data[name]:
                mpl_smin = 18
                mpl_smax = 180
                kwargs['s'] = (mpl_smax - mpl_smin) * self.data[name]['size_scale'] + mpl_smin
            else:
                kwargs['s'] = Visualization.MPL_SIZE
        return kwargs

    def plot_scatter(self, plot_args=None, scatter_args=None, multiple_graphs=False, title=None, limits=None):
        """
        Visualize data using maplotlib scatter plots.

        :param plot_args: additional keyword arguments passed to :py:class:`matplotlib.pyplot.subplots`
        :type plot_args: dict

        :param scatter_args: dict of additional keyword arguments passed to \
        :py:class:`matplotlib.pyplot.scatter` for each group.
        :type scatter_args: dist

        :param multiple_graphs: if each data entry should be ploted separately or on the same plot. \
        If ploted separately a shape may be provided.
        :type multiple_graphs: bool | tuple(nrows, ncols)
        """
        if mpl_import_err is not None:
            raise VisualizationError("Missing MatPlotLib dependency.") from mpl_import_err

        plot_args = _validate_mplplot_args(plot_args)
        if multiple_graphs:
            plot_args['nrows'], plot_args['ncols'] = self.__get_grid_shape(multiple_graphs)

        fig, axs = plt.subplots(**plot_args)
        if multiple_graphs:
            axs = axs.flatten()
            if len(self.data.keys())%plot_args['ncols'] != 0:
                fig.delaxes(axs[-1])
                axs = numpy.delete(axs, -1)
            data_keys = list(self.data.keys())
            for index, ax in enumerate(axs):
                name = data_keys[index]
                coords, axes_labels = _get_coords(self.data[name]["points"])

                group_args = self.__validate_mplscatter_args(scatter_args, name)
                if "cmap" in group_args:
                    group_args['c'] = self.data[name]["data"]

                ax.scatter(*coords, label=name, **group_args)
                ax.set_xlabel(axes_labels[0])
                ax.set_xlim(limits[0])
                ax.set_ylabel(axes_labels[1])
                ax.set_ylim(limits[1])
                ax.grid(linestyle='--', linewidth=1)
                ax.legend(loc='upper right', fontsize=12)
                if title is not None:
                    ax.set_title(title)
        else:
            for index, (name, data) in enumerate(self.data.items()):
                x_coords = list(map(lambda point: point[0], data["points"]))
                y_coords = list(map(lambda point: point[1], data["points"]))

                group_args = self.__validate_mplscatter_args(scatter_args, name)
                if "cmap" in group_args:
                    group_args['c'] = data["data"]

                axs.scatter(x_coords, y_coords, label=name, **group_args)
                axs.grid(linestyle='--', linewidth=1)
                axs.legend(loc='upper right', fontsize=12)
                if title is not None:
                    axs.set_title(title)

            plt.axis('scaled')
