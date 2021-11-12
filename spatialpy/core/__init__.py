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

import logging
from .BoundaryCondition import *
from .cleanup import *
from .DataFunction import *
from .Domain import *
from .Geometry import *
from .InitialCondition import *
from .Model import *
from .Parameter import *
from .Reaction import *
from .Result import *
from .spatialpyError import *
from .Species import *
from .VTKReader import *
from spatialpy.__version__ import __version__

_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
version = __version__
log = logging.getLogger()
log.setLevel(logging.WARNING)
log.addHandler(_handler)

__all__ = [s for s in dir() if not s.startswith('_')]