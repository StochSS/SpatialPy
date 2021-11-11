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

# Base Module Expections
class ModelError(Exception):
    pass


class ResultError(Exception):
    pass


class VTKReaderError(Exception):
    """Base class for exceptions in VTKReader module."""
    pass


# Model Component Exceptions
class BoundaryConditionError(ModelError):
    pass


class DataFunctionError(ModelError):
    pass


class DomainError(ModelError):
    pass


class GeometryError(ModelError):
    pass


class InitialConditionError(ModelError):
    pass


class ParameterError(ModelError):
    pass


class ReactionError(ModelError):
    pass


class SpeciesError(ModelError):
    pass


# Result Exceptions


# VTKReader Exceptions
class VTKReaderIOError(VTKReaderError):
    """Exception raised for I/O errors."""
    def __init__(self, message):
        self.message = message