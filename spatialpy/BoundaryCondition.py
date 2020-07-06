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
                 type_id=None,
                 species=None,
                 deterministic=True,
                 property=None,
                 value=None):
        """
            Conditions (one or more of the following must be set):
                 xmin, xmax: (float) min or max value in the x dimension
                 ymin, ymax: (float) min or max value in the y dimension
                 zmin, zmax: (float) min or max value in the z dimension
                 type_id: type (subdomain) of the partciles
            Targets (one of the following must be set):
                property: (str), 'nu', 'rho','v'
                spesicies: (str) name of a chemical species.  
                           Must also set deterministic=True/False flag.
            Assignment:
                value: (float or float[3]), value property will take in region 
                       defined by the conditions

        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.type_id = type_id
        self.species = None
        self.property = None
        self.deterministic = deterministic

        if( species is not None and property is not None):
            raise ModelError("Can not set both species and property")
        if species is not None:
            self.species = species
            self.deterministic = deterministic
        else:
            self.property = property

        if value is None:
            raise ModelError("Must set value")
        self.value = value

        


