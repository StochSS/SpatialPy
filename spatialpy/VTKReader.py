import numpy


class VTKReader:
    """VTKReader.py: SpatialPy minimal VTK legacy file reader."""
    """Reference: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf"""

    def __init__(self):

        self.filename = None
        self.pointdatatype = None
        self.numpoints = None
        self.points = None
        self.arrays = None
        self.datatypes = {
            "int": "int32",
            "float": "float32",
            "double": "float64",
        }

    def setfilename(self, filename):
        """Set filename.
        Args:
            (str) filename
        """

        self.filename = filename

    def getnumarrays(self):
        """Get (int) number of arrays."""

        length = len(self.arrays)

        if length > 0:
            return length
        else:
            return None

    def getarrayname(self, i):
        """Get (str) array name.
        Args:
            (int) index"""

        arrayids = list(self.arrays.keys())

        if i <= len(arrayids):
            return arrayids[i]
        else:
            return None

    def getarrays(self):
        """Get (dict) of arrays."""

        return self.arrays

    def getnumpoints(self):
        """Get (int) number of points."""

        return self.numpoints

    def getpoints(self):
        """Get (list) points."""
        return self.points

    def readnumeric(self, fd):
        """Read numeric data to be returned in a list.
        Args:
            file descriptor
        Return:
            (list) numeric data read"""

        numericlist = []

        for line in fd:
            if line.lstrip("-")[:1].isnumeric():
                numericlist.extend(line.strip().split())
            else:
                break

        return numericlist

    def readpoints(self, fd, row, col, datatype):
        """Read numeric data to be returned in a shaped numpy array.
        Args:
            file descriptor
            (int) number of rows
            (int) number of columns
            (str) VTK datatype
        Return:
            (numpy array) of shape row by col of (numpy) dtype
        """

        pointlist = numpy.array(self.readnumeric(fd), dtype=self.datatypes[datatype])
        if col > 1:
            pointlist = pointlist.reshape(row, col)

        return pointlist

    # This should be refactored at some point and probably eliminated
    # In factor of readnumeric()
    def populatearrays(self, vtkdata, arraydata, col, row, name, datatype):
        """Populate arrays with data to be added to dict of shaped arrays.
        Args:
            (dict) array data
            (array) array data to be added to dict
            (int) number of rows
            (int) number of columns
            (str) name of dict key
            (str) VTK datatype
        """

        array = numpy.array(arraydata, dtype=datatype)
        if col > 1:
            array = array.reshape(row, col)
        vtkdata[name] = array

    def readarrays(self, fd):
        """Read all arrays and data in file.
        Args:
            file descriptor
        Return:
            (dict) array data
        """
        vtkdata = {}
        arraydata = []

        name, col, row, datatype = fd.readline().strip().split()
        col = int(col)
        row = int(row)

        for line in fd:

            # Array names MUST begin with a letter
            if line[:1].isalpha():

                self.populatearrays(vtkdata, arraydata, col, row, name, datatype)
                arraydata = []

                name, col, row, datatype = line.strip().split()
                col = int(col)
                row = int(row)

            elif line.lstrip("-")[:1].isnumeric():
                arraydata.extend(line.strip().split())

            elif line.isspace():
                continue

        self.populatearrays(vtkdata, arraydata, col, row, name, datatype)

        return vtkdata

    def readfile(self):
        """Read VTK file."""

        with open(self.filename) as fd:
            fd.readline()
            fd.readline()

            # We only output ASCII so we can ignore BINARY
            if fd.readline().strip().upper() != "ASCII":
                raise VTKReaderIOError("{0} doesn't look like a valid ASCII VTK file.".format(self.filename))

            fd.readline()

            _, self.numpoints, self.pointdatatype = fd.readline().strip().split()
            self.numpoints = int(self.numpoints)

            self.points = self.readpoints(fd, self.numpoints, 3, self.pointdatatype)

            for line in fd:
                if line[:5] != "FIELD":
                    continue
                else:
                    break

            self.arrays = self.readarrays(fd)


class VTKReaderError(Exception):
    """Base class for exceptions in VTKReader module."""

    pass


class VTKReaderIOError(VTKReaderError):
    """Exception raised for I/O errors."""

    def __init__(self, message):
        self.message = message
