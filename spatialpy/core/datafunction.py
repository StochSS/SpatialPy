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
from spatialpy.core.spatialpyerror import DataFunctionError

class DataFunction():
    """
    Abstract class used to constuct the data function.

    :param name: Name of the Data Function.
    :type name: str

    :raises DataFunctionError: If a name is not provided.
    """

    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise DataFunctionError("DataFunction must have a 'name'")

    def map(self, point):
        """
        This method must be overridden by the DataFunction subclass.

        NOTE: The spatial location is evaulated at t=0 and is not \
              re-evaluated as the fluid domain moves over time.

        :param point: The x, y, z position.
        :type point: float[3]

        :returns: Value of function at this spatial location.
        :rtype: float
        """
        raise DataFunctionError(f"{self.name}: DataFunction.map() must be implemented.")
