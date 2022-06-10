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
import numpy as np
from collections.abc import Iterator

from .spatialpyerror import TimespanError

class TimeSpan(Iterator):
    """
    Model timespan that describes the duration to run the simulation and at which timepoint to sample
    the species populations during the simulation.

    :param items: Evenly-spaced list of times at which to sample the species populations during the simulation. 
            Best to use the form np.linspace(<start time>, <end time>, <number of time-points, inclusive>)
    :type items: list, tuple, range, or numpy.ndarray
    
    :param timestep_size: Size of each timestep in seconds
    :type timestep_size: int | float

    :raises TimespanError: items is an invalid type.
    """
    def __init__(self, items, timestep_size=None):
        if isinstance(items, (list, tuple, range)):
            items = np.array(items)

        self.validate(items=items, timestep_size=timestep_size)

        if timestep_size is None:
            timestep_size = float(items[1] - items[0])
        self.timestep_size = timestep_size

        items_diff = np.diff(items)
        self._set_timesteps(items_diff[0], len(items_diff))

        self.validate(coverage="initialized")

    def __eq__(self, o):
        return self.items.__eq__(o).all()

    def __getitem__(self, key):
        return self.items.__getitem__(key)

    def __iter__(self):
        return self.items.__iter__()

    def __len__(self):
        return self.items.__len__()

    def __next__(self):
        return self.items.__next__()

    def __str__(self):
        return self.items.__str__()

    def _ipython_display_(self):
        print(self)

    def _set_timesteps(self, output_interval, num_steps):
        if self.timestep_size is None:
            self.timestep_size = output_interval

        self.output_freq = output_interval / self.timestep_size
        self.num_timesteps = math.ceil(num_steps * self.output_freq)

        output_steps = np.arange(
            0, self.num_timesteps + self.timestep_size, self.output_freq
        )
        self.output_steps = np.unique(np.round(output_steps).astype(int))
        self.items = np.zeros((self.output_steps.size), dtype=float)
        for i, step in enumerate(self.output_steps):
            self.items[i] = step * self.timestep_size

    @classmethod
    def linspace(cls, t=20, num_points=None, timestep_size=None):
        """
        Creates a timespan using the form np.linspace(0, <t>, <num_points, inclusive>).

        :param t: End time for the simulation.
        :type t: float | int

        :param num_points: Number of sample points for the species populations during the simulation.
        :type num_points: int

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: int | float

        :returns: Timespan for the model.
        :rtype: spatialpy.TimeSpan

        :raises TimespanError: t or num_points are None, <= 0, or invalid type.
        """
        if t is None or not isinstance(t, (int, float)) or t <= 0:
            raise TimespanError("t must be a positive float or int.")
        if num_points is not None and (not isinstance(num_points, int) or num_points <= 0):
            raise TimespanError("num_points must be a positive int.")

        if num_points is None:
            num_points = int(t / 0.05) + 1
        items = np.linspace(0, t, num_points)
        return cls(items, timestep_size=timestep_size)

    @classmethod
    def arange(cls, increment, t=20, timestep_size=None):
        """
        Creates a timespan using the form np.arange(0, <t, inclusive>, <increment>).

        :param increment: Distance between sample points for the species populations during the simulation.
        :type increment: float | int

        :param t: End time for the simulation.
        :type t: float | int

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: int | float

        :returns: Timespan for the model.
        :rtype: spatialpy.TimeSpan

        :raises TimespanError: t or increment are None, <= 0, or invalid type.
        """
        if t is None or not isinstance(t, (int, float)) or t <= 0:
            raise TimespanError("t must be a positive floar or int.")
        if not isinstance(increment, (float, int)) or increment <= 0:
            raise TimespanError("increment must be a positive float or int.")

        items = np.arange(0, t + increment, increment)
        return cls(items, timestep_size=timestep_size)

    def validate(self, items=None, timestep_size=None, coverage="build"):
        """
        Validate the models time span

        :param timestep_size: Size of each timestep in seconds
        :type timestep_size: int | float

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises TimespanError: Timespan is an invalid type, empty, not uniform, contains a single \
                               repeated value, or contains a negative initial time.
        """
        if coverage in ("all", "build"):
            if hasattr(self, "items") and items is None:
                items = self.items

            if not isinstance(items, np.ndarray):
                if not isinstance(items, (list, tuple, range)):
                    raise TimespanError("Timespan must be of type: list, tuple, range, or numpy.ndarray.")
                items = np.array(items)
                if items is not None:
                    self.items = items

            if len(items) == 0:
                raise TimespanError("Timespans must contain values.")
            if items[0] < 0:
                raise TimespanError("Simulation must run from t=0 to end time (t must always be positive).")
            
            first_diff = items[1] - items[0]
            other_diff = items[2:] - items[1:-1]
            isuniform = np.isclose(other_diff, first_diff).all()

            if coverage == "build" and not isuniform:
                raise TimespanError("SpatialPy only supports uniform timespans.")
            if first_diff == 0 or np.count_nonzero(other_diff) != len(other_diff):
                raise TimespanError("Timespan can't contain a single repeating value.")

        if coverage in ("all", "build", "timestep_size"):
            if hasattr(self, "timestep_size") and timestep_size is None:
                timestep_size = self.timestep_size

            if timestep_size is not None:
                if not isinstance(timestep_size, (int, float)):
                    raise TimespanError("timestep_size must be of type int or float.")
                if timestep_size <= 0:
                    raise TimespanError("timestep_size must be a positive value.")

        if coverage in ("all", "initialized"):
            if self.timestep_size is None:
                raise TimespanError("timestep_size can't be None type.")
            if self.output_freq is None:
                raise TimespanError("output_freq can't be None type.")
            if not isinstance(self.output_freq, (int, float)):
                raise TimespanError("output_freq must be of type int or float.")
            if self.output_freq < self.timestep_size:
                raise TimespanError("timestep_size exceeds output_frequency.")
            if self.num_timesteps is None:
                raise TimespanError("num_timesteps can't be None type.")
            if not isinstance(self.num_timesteps, int):
                raise TimespanError("num_timesteps must be of type int.")
            if self.num_timesteps <= 0:
                raise TimespanError("num_timesteps must be a positive int.")
            if self.output_steps is None:
                raise TimespanError("output_steps can't be None type.")
            if not isinstance(self.output_steps, (np.ndarray)):
                raise TimespanError("output_steps must be of type numpy.ndarray.")
            if self.items.size != self.output_steps.size:
                raise TimespanError("output_steps must be the same size as items.")
