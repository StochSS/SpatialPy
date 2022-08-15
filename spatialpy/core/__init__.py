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

import logging
from spatialpy.__version__ import __version__
from .boundarycondition import BoundaryCondition
from .cleanup import (
	cleanup_tempfiles,
	cleanup_core_files,
	cleanup_build_files,
	cleanup_result_files
)
from .datafunction import DataFunction
from .domain import Domain
from .geometry import (
	CombinatoryGeometry,
	Geometry,
	GeometryAll,
	GeometryExterior,
	GeometryInterior
)
from .initialcondition import (
	InitialCondition,
	PlaceInitialCondition,
	UniformInitialCondition,
	ScatterInitialCondition
)
from .lattice import (
	CartesianLattice,
	SphericalLattice,
	CylindricalLattice,
	XMLMeshLattice,
	MeshIOLattice,
	StochSSLattice
)
from .model import Model, export_StochSS
from .parameter import Parameter
from .reaction import Reaction
from .result import Result
from .spatialpyerror import *
from .species import Species
from .timespan import TimeSpan
from .transformation import (
	Transformation,
	TranslationTransformation,
	RotationTransformation,
	ReflectionTransformation,
	ScalingTransformation
)
from .visualization import Visualization
from .vtkreader import *

_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
version = __version__
log = logging.getLogger()
log.setLevel(logging.WARNING)
log.addHandler(_handler)

__all__ = [s for s in dir() if not s.startswith('_')]
