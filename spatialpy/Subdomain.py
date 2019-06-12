

class SubDomain:
    """ SubDomain class provides a method for taging parts of the spatial domain as seperate parts"""

    def __init__(self):
        pass


    def inside(self, x, on_boundary):
        raise Exception("Subclasses of spatialpy.SubDomain must implement the inside() metho")



