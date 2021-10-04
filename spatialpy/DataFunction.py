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


class DataFunction():
    """ Abstract class used to constuct the data function. """

    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise Exception("DataFunction must have a 'name'")

    def map(self, x):
        """
        This method must be overridden by the DataFunction subclass.

        NOTE: The spatial location is evaulated at t=0 and is not 
              reevaluated as the fluid domain moves over time.

        Args:
            x (vector of 3 doubles): the x,y,z position to evalute the function.
        Returns:
            float: value of the data function at the spatial location given by 'x'.
        """
        raise Exception("DataFunction.map() must be implemented.")


