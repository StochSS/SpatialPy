#__all__=['model','spatialpy']
from .spatialpy import *

import sys
if (sys.version_info < (3,0)):
    raise Exception("SpatialPy only works in Python 3.0 and higher")


#from spatialpy.core import *
#import spatialpy.sbml
#import spatialpy.solvers
#import spatialpy.example_models
