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
import numpy

from spatialpy.core.spatialpyerror import InitialConditionError

class InitialCondition():
    """
    Class used to defined initial conditions in SpatialPy.
    SubClasses must implement the 'apply(model)' method, which
    direction modifies the model.u0[species, voxel] matrix.
    """
    def apply(self, model):
        """
        Set the initial condition of the species to the count.

        :param model: Model contianing the target species.
        :type model: spatialpy.core.model.Model
        """
        raise InitialConditionError("spatialpy.InitialCondition subclasses must implement apply()")

class PlaceInitialCondition(InitialCondition):
    """
    Class used to defined the place initial condition in SpatialPy.

    :param species: The species to set the initial condition.
    :type species: spatialpy.core.species.Species

    :param count: The initial condition for the target species.
    :type count: int

    :param location: X, Y, Z coordinates to place the initial condition.
    :type location: float[3]
    """
    def __init__(self, species, count, location):
        self.species = species
        self.count = count
        self.location = location

    def __str__(self):
        print_string = f"{self.species.name}: {str(self.count)}, at: {str(self.location)}"
        return print_string

    def apply(self, model):
        """
        Set the initial condition of the species to the count at the location.

        :param model: Model contianing the target species.
        :type model: spatialpy.core.model.Model
        """
        if isinstance(self.species, str):
            if self.species not in model.listOfSpecies:
                raise InitialConditionError(f"Species {self.species} does not exist in the model.")
            self.species = model.listOfSpecies[self.species]
        spec_name = self.species.name
        spec_ndx = None
        for index, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species:
                spec_ndx = index
                break
        vtx = model.domain.closest_vertex(self.location)
        model.u0[spec_ndx, vtx] += self.count

class UniformInitialCondition(InitialCondition):
    """
    Class used to defined the uniform initial condition in SpatialPy.

    :param species: The species to set the initial condition.
    :type species: spatialpy.core.species.Species

    :param count: The initial condition for the target species.
    :type count: int

    :param types: Types of the particles to place the initial condition.
    :type types: list
    """
    def __init__(self, species, count, types=None):
        self.species = species
        self.count = count
        if types is None:
            self.types = types
        else:
            self.types = []
            for type_id in types:
                self.types.append(f"type_{type_id}")

    def __str__(self):
        print_string = f"{self.species.name}: {str(self.count)}, Uniformly distrbuted in: {str(self.types)}"
        return print_string

    def apply(self, model):
        """
        Set 'count' of 'species' in over the list of types (all types if None).

        :param model: Model contianing the target species.
        :type model: spatialpy.core.model.Model
        """
        if isinstance(self.species, str):
            if self.species not in model.listOfSpecies:
                raise InitialConditionError(f"Species {self.species} does not exist in the model.")
            self.species = model.listOfSpecies[self.species]
        spec_name = self.species.name
        spec_ndx = None
        for index, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species:
                spec_ndx = index
                break

        if self.types is None:
            nvox = model.domain.get_num_voxels()
            for vtx in range(nvox):
                model.u0[spec_ndx, vtx] += self.count
        else:
            for i in range(model.domain.get_num_voxels()):
                type_id = model.domain.type_id[i]
                if type_id in self.types:
                    model.u0[spec_ndx, i] += self.count


class ScatterInitialCondition(InitialCondition):
    """
    Class used to defined the scatter initial condition in SpatialPy.

    :param species: The species to set the initial condition.
    :type species: spatialpy.core.species.Species

    :param count: The initial condition for the target species.
    :type count: int

    :param types: Types of the particles to place the initial condition.
    :type types: list
    """
    def __init__(self, species, count, types=None):
        self.species = species
        self.count = count
        if types is None:
            self.types = types
        else:
            self.types = []
            for type_id in types:
                self.types.append(f"type_{type_id}")

    def __str__(self):
        print_string = f"{self.species.name}: {str(self.count)}, Scatter over: {str(self.types)}"
        return print_string

    def apply(self, model):
        """
        Scatter 'count' of 'species' randomly over the list of types (all types if None).

        :param model: Model contianing the target species.
        :type model: spatialpy.core.model.Model
        """
        if isinstance(self.species, str):
            if self.species not in model.listOfSpecies:
                raise InitialConditionError(f"Species {self.species} does not exist in the model.")
            self.species = model.listOfSpecies[self.species]
        spec_name = self.species.name
        spec_ndx = None
        for index, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species:
                spec_ndx = index
                break

        if self.types is None:
            nvox = model.domain.get_num_voxels()
            for _ in range(self.count):
                vtx = numpy.random.randint(0, nvox)
                model.u0[spec_ndx, vtx] += 1
        else:
            allowed_voxels = []
            for i in range(model.domain.get_num_voxels()):
                type_id = model.domain.type_id[i]
                if type_id in self.types:
                    allowed_voxels.append(i)
            nvox = len(allowed_voxels)
            if nvox==0:
                message = "ScatterInitialCondition has zero voxels to scatter in. "
                message += f"Species={self.species.name} count={self.count} types={self.types}"
                raise InitialConditionError(message)
            for _ in range(self.count):
                v_ndx = numpy.random.randint(0, nvox)
                vtx = allowed_voxels[v_ndx]
                model.u0[spec_ndx, vtx] += 1
