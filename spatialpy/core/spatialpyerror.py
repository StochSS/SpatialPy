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

# Base Module Expections
class ModelError(Exception):
    """
    Class for exceptions in the model module.
    """

class ResultError(Exception):
    """
    Class for exceptions in the results module.
    """

class VisualizationError(Exception):
    """
    Class for exceptions in the visualization module.
    """

class VTKReaderError(Exception):
    """
    Bass class for exceptions in the vtkreader module.
    """

class SimulationError(Exception):
    """
    Class for exceptions in the simulation module.
    """

# Model Component Exceptions
class BoundaryConditionError(ModelError):
    """
    Base class for exceptions in the boundarycondition module.
    """

class DataFunctionError(ModelError):
    """
    Class for exceptions in the datafunction module.
    """

class DomainError(ModelError):
    """
    Class for exceptions in the domain module.
    """

class GeometryError(ModelError):
    """
    Class for exceptions in the geometry module.
    """

class InitialConditionError(ModelError):
    """
    Class for exceptions in initailcondition module.
    """

class LatticeError(ModelError):
    """
    Class for exceptions in lattice module.
    """

class ParameterError(ModelError):
    """
    Class for exceptions in parameter module.
    """

class ReactionError(ModelError):
    """
    Class for exceptions in reaction module.
    """

class SpeciesError(ModelError):
    """
    Class for exceptions in the species module.
    """

class TimespanError(ModelError):
    """
    Class for exceptions in the timespan module.
    """

class TransformationError(ModelError):
    """
    Class for exceptions in the transformation module.
    """

# Result Exceptions


# Visualization Exceptions


# VTKReader Exceptions
class VTKReaderIOError(VTKReaderError):
    """
    Exception raised for I/O errors.
    """
    def __init__(self, message):
        super().__init__(message)
        self.message = message

# Simulation Exceptions
class SimulationTimeout(SimulationError):
    """
    Exception raised for timeout errors.
    """
