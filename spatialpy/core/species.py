'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2019 - 2022 SpatialPy developers.

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
from spatialpy.core.parameter import Parameter
from spatialpy.core.spatialpyerror import SpeciesError

class Species():
    """
    Model of a biochemical species. Must be assigned a diffusion coefficent.

    :param name: Name of the Species
    :type name: str

    :param diffusion_coefficient: Non-constant coefficient of diffusion for Species.
    :type diffusion_coefficient: float

    :param restrict_to: Set the diffusion coefficient to zero for 'species' in all types not in 'listOfTypes'.
    This effectively restricts the movement of 'species' to the types specified in 'listOfTypes'.
    :type restrict_to: int, str, list of ints or list of strs
    """
    def __init__(self, name=None, diffusion_coefficient=None, restrict_to=None):
        if name is None:
            raise SpeciesError("Species must have a name")
        if not isinstance(name, str):
            raise SpeciesError("Species name must be a string")

        if  diffusion_coefficient is None:
            raise SpeciesError("Species must have a diffusion_coefficient.")
        if not (isinstance(diffusion_coefficient, (Parameter, float, int)) or \
                    type(diffusion_coefficient).__name__ == 'Parameter'):
            raise SpeciesError("Diffusion coefficient must be a spatialpy.Parameter, float, or int.")
        if isinstance(diffusion_coefficient, (float, int)) and diffusion_coefficient < 0:
            raise SpeciesError("Diffusion coefficient must be non-negative.")

        if not (restrict_to is None or isinstance(restrict_to, (str, int, list))):
            raise SpeciesError("Restrict_to must be an int, str or list of ints or strs.")
        if restrict_to is not None and isinstance(restrict_to, (int, str)):
            restrict_to = [restrict_to]

        self.name = name
        self.diffusion_coefficient = diffusion_coefficient
        if restrict_to is None:
            self.restrict_to = restrict_to
        else:
            self.restrict_to = []
            for type_id in restrict_to:
                self.restrict_to.append(f"type_{type_id}")

    def __str__(self):
        print_string = f"{self.name}: {str(self.diffusion_coefficient)}"
        return print_string

    def set_diffusion_coefficient(self, diffusion_coefficient):
        """
        Setter method for non-constant coefficient of diffusion for Species.

        :param diffusion_coefficient: Non-constant coefficient of diffusion for Species.
        :type diffusion_coefficient: float

        :raises SpeciesError: If diffusion_coefficient is negative.
        """
        if not (isinstance(diffusion_coefficient, (Parameter, float, int)) or \
                            type(diffusion_coefficient).__name__ == 'Parameter'):
            raise SpeciesError("Diffusion coefficient must be a spatialpy.Parameter, float, or int.")
        if diffusion_coefficient < 0:
            raise SpeciesError("Diffusion coefficient must be non-negative.")

        self.diffusion_coefficient = diffusion_coefficient
