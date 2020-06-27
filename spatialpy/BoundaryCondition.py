import numpy
from spatialpy.Model import ModelError


class BoundaryCondition():
    """ Set spatial regions of the domain where a property of 
        particles are held constant (updated each simulation step)
    """
    def __init__(self, 
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 zmin=None, zmax=None,
                 species=None,
                 property=None,
                 value=None):
        """
            Arguments:

        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.species =  species
        self.property = property
        self.value = value
        if( species is not None and property is not None):
            ModelError("Can not set both species and property")
        if value is None:
            ModelError("Must set value")


