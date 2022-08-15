# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2022 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy

from spatialpy.core.spatialpyerror import VTKReaderIOError

def _is_valid_num(numstr):
    try:
        float(numstr)
    except ValueError:
        return False

    return True

def _read_arrays(data_file):
    vtkdata = {}
    arraydata = []

    for line in data_file:
        if line.isspace():
            continue
        try:
            name, col, row, datatype = line.strip().split()
            # got another PointData section
        except Exception as err:
            print(f"Error: {err}")
            print("on line >>>")
            print(line, end='')
            print(f"<<< {data_file.tell()}")
            raise err
        col = int(col)
        row = int(row)

        # now read row*col number of values
        arraydata = []
        while len(arraydata) < row * col:
            line = data_file.readline()
            arraydata.extend(line.strip().split())

        _populate_arrays(vtkdata, arraydata, col, row, name, datatype)

    return vtkdata

def _read_numeric(data_file):
    numericlist = []

    for line in data_file:
        parts = line.strip().split()
        if len(parts) == 0:
            break
        if _is_valid_num(parts[0]):
            numericlist.extend(parts)
        else:
            break

    return numericlist

# This should be refactored at some point and probably eliminated
# In factor of __read_numeric()
def _populate_arrays(vtkdata, arraydata, col, row, name, datatype):
    array = numpy.array(arraydata, dtype=datatype)
    if col > 1:
        array = array.reshape(row, col)
    vtkdata[name] = array

class VTKReader:
    """
    VTKReader.py: SpatialPy minimal VTK legacy file reader.
    Reference: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

    :param filename: name of SpatialPy VTK file
    :type filename: str

    :param debug: If true, will print debugging information.
    :type debug: bool
    """
    def __init__(self, filename=None, debug=False):

        self.filename = filename
        self.pointdatatype = None
        self.numpoints = None
        self.points = None
        self.arrays = None
        self.debug = debug
        self.datatypes = {
            "int": "int32",
            "float": "float32",
            "double": "float64",
        }

    def __read_points(self, data_file, row, col, datatype):
        pointlist = numpy.array(_read_numeric(data_file), dtype=self.datatypes[datatype])
        if col > 1:
            pointlist = pointlist.reshape(row, col)

        return pointlist

    def set_filename(self, filename):
        """
        Set the filename.

        :params filename: Filename
        :type filename: str
        """
        self.filename = filename

    def get_array_name(self, index):
        """
        Get the array name.

        :param index: index
        :type index: int

        :returns: Name/Key for arraid id
        :rtype: int | None
        """
        arrayids = list(self.arrays.keys())

        if index <= len(arrayids):
            return arrayids[index]
        return None

    def get_arrays(self):
        """
        Get the dictionary of arrays.

        :returns: dictionary of data arrays
        :rtype: dict
        """
        return self.arrays

    def get_num_points(self):
        """
        Get the number of points.

        :returns: Number of points
        :rtype: int
        """
        return self.numpoints

    def get_points(self):
        """
        Get the list of points.

        :returns: List of points.
        :rtype: list
        """
        return self.points

    def read_file(self):
        """
        Read VTK file.

        :raises VTKReaderIOError: Invalid ASCII VTK file
        """
        with open(self.filename, encoding="utf-8") as data_file:
            if self.debug:
                print(f"open({self.filename})")
            tmp =  data_file.readline()
            tmp =  data_file.readline()

            # We only output ASCII so we can ignore BINARY
            tmp =  data_file.readline()
            if tmp.strip().upper() != "ASCII":
                raise VTKReaderIOError(f"{self.filename} doesn't look like a valid ASCII VTK file.")

            tmp =  data_file.readline()

            tmp =  data_file.readline()
            _, self.numpoints, self.pointdatatype = tmp.strip().split()
            self.numpoints = int(self.numpoints)

            self.points = self.__read_points(data_file, self.numpoints, 3, self.pointdatatype)

            for line in data_file:
                if line[:5] != "FIELD":
                    continue
                break

            self.arrays = _read_arrays(data_file)
