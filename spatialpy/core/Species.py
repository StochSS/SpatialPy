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
from spatialpy.core.spatialpyError import ModelError

class Species():
    """ Model of a biochemical species. Must be assigned a diffusion coefficent.

            :param name: Name of the Species
            :type name: str
            :param diffusion_coefficient: non-constant coefficient of diffusion for Species
            :type diffusion_coefficient: float
            """

    reserved_names = ["x", "vol","sd","data_fn","t","debug_flag","Spatialpy"]



    def __init__(self,name=None, diffusion_coefficient=None):
        # A species has a name (string) and an initial value (positive integer)
        if name is None:
            raise ModelError("Species must have a name")
        else:
            self.name = name
        if  diffusion_coefficient is not None:
            self.diffusion_coefficient=diffusion_coefficient
        else:
            raise ModelError("Species must have a diffusion_coefficient")


    def __str__(self):
        print_string = f"{self.name}: {str(self.diffusion_coefficient)}"
        return print_string
