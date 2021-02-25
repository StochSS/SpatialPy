class Geometry:
    """ Geometry class provides a method for tagging parts of the spatial domain as separate parts"""

    def __init__(self):
        pass

    def inside(self, x, on_boundary):
        raise Exception("Subclasses of spatialpy.Geometry must implement the inside() method")
