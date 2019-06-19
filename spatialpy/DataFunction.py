import spatialpy


class DataFunction():
    """ Abstract class used to constuct the data vector. """
    name = None
    def __init__(self, name=None):
        if name is not None:
            self.name = name
        if self.name is None:
            raise Exception("DataFunction must have a 'name'")

    def map(self, x):
        """ map() takes the coordinate 'x' and returns a double.
        Args:
            x: a list of 3 ints.
        Returns:
            A float
        """
        raise Exception("DataFunction.map() not implemented.")


