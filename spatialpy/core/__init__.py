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

import logging
from spatialpy.__version__ import __version__
from .boundarycondition import *
from .cleanup import *
from .datafunction import *
from .domain import *
from .geometry import *
from .initialcondition import *
from .model import Model
from .parameter import *
from .reaction import *
from .result import *
from .spatialpyerror import *
from .species import *
from .vtkreader import *

_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
version = __version__
log = logging.getLogger()
log.setLevel(logging.WARNING)
log.addHandler(_handler)

__all__ = [s for s in dir() if not s.startswith('_')]
