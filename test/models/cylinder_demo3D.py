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

import os
import numpy
import sys
from collections import OrderedDict
sys.path.append("../..")
import spatialpy
# Global Constants
MAX_X_DIM = 5.0
MIN_X_DIM = -5.0
TOL = 1e-9

class Edge1(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return abs(x[0] - MAX_X_DIM) < 0.05
class Edge2(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return abs(x[0] - MIN_X_DIM) < 0.05
class Middle(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return abs(x[0] - MIN_X_DIM) >= 0.05
class cylinderDemo3D(spatialpy.Model):
    EDGE1 = Edge1.__name__
    EDGE2 = Edge2.__name__
    MIDDLE = Middle.__name__

    def __init__(self, model_name="cylinder_demo3d"):
        spatialpy.Model.__init__(self, model_name)

        self.timestep_size = 1

        # System constants
        D_const = 0.1

        # Define Species
        A = spatialpy.Species(name="A", diffusion_coefficient=D_const, restrict_to=[self.MIDDLE, self.EDGE1])
        B = spatialpy.Species(name="B", diffusion_coefficient=D_const, restrict_to=[self.MIDDLE, self.EDGE2])
        self.add_species([A, B])

        # Define Geometry
        # Make sure that we have the correct path to the mesh file even if we are not executing from the basedir.
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.domain = spatialpy.Domain.read_xml_mesh(basedir + 'cylinder.xml')

        # Define Types
        self.set_properties(Middle(), self.MIDDLE)
        self.set_properties(Edge1(), self.EDGE1)
        self.set_properties(Edge2(), self.EDGE2)

        vol = self.domain.get_vol()
        type_id = self.domain.type
        left = numpy.sum(vol[type_id == self.domain.get_type_def(self.EDGE1)])
        right = numpy.sum(vol[type_id == self.domain.get_type_def(self.EDGE2)])
        print("left "+str(left)+" right "+str(right))

        k_react = spatialpy.Parameter(name="k_react", expression=1.0)
        k_creat1 = spatialpy.Parameter(name="k_creat1",
                                     expression=100/left)
        k_creat2 = spatialpy.Parameter(name="k_creat2",
                                     expression=100/right)
        self.add_parameter([k_react, k_creat1, k_creat2])


        # Define Reactions
        R1 = spatialpy.Reaction(reactants=None, products={A:1},
                                rate=k_creat1, restrict_to=self.EDGE1)
        R2 = spatialpy.Reaction(reactants=None, products={B:1},
                              rate=k_creat2, restrict_to=self.EDGE2)
        R3 = spatialpy.Reaction(reactants={A:1, B:1}, products=None,
                              rate=k_react)
        self.add_reaction([R1, R2, R3])

        # Define simulation timespan
        #self.set_timesteps(1, 200)
        self.timespan(range(500))
