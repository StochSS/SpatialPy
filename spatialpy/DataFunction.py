import spatialpy



class DataFunction(spatialpy.BoundaryCondition):
    """ Abstract class used to constuct the data function. """

    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise Exception("DataFunction must have a 'name'")

    def expression(self):
        """ 
        """
        raise Exception("DataFunction.expression() must be implemented.")


