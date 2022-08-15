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

def create_cylinder_demo_3D(model_name="cylinder_demo3d", parameter_values=None):
    model = spatialpy.Model(model_name)

    # Set domain type ids
    model.EDGE1 = Edge1.__name__
    model.EDGE2 = Edge2.__name__
    model.MIDDLE = Middle.__name__

    # System constants
    D_const = 0.1

    # Define Geometry
    # Make sure that we have the correct path to the mesh file even if we are not executing from the basedir.
    basedir = os.path.dirname(os.path.abspath(__file__))
    domain = spatialpy.Domain.read_xml_mesh(f'{basedir}cylinder.xml')
    
    # Define Types
    domain.set_properties(Middle(), model.MIDDLE)
    domain.set_properties(Edge1(), model.EDGE1)
    domain.set_properties(Edge2(), model.EDGE2)

    model.add_domain(domain)

    # Define Species
    A = spatialpy.Species(name="A", diffusion_coefficient=D_const, restrict_to=[model.MIDDLE, model.EDGE1])
    B = spatialpy.Species(name="B", diffusion_coefficient=D_const, restrict_to=[model.MIDDLE, model.EDGE2])
    model.add_species([A, B])

    vol = model.domain.get_vol()
    type_id = model.domain.type
    left = numpy.sum(vol[type_id == model.domain.get_type_def(model.EDGE1)])
    right = numpy.sum(vol[type_id == model.domain.get_type_def(model.EDGE2)])
    print(f"left {left} right {right}")

    k_react = spatialpy.Parameter(name="k_react", expression=1.0)
    k_creat1 = spatialpy.Parameter(name="k_creat1", expression=100 / left)
    k_creat2 = spatialpy.Parameter(name="k_creat2", expression=100 / right)
    model.add_parameter([k_react, k_creat1, k_creat2])

    # Define Reactions
    R1 = spatialpy.Reaction(reactants=None, products={A: 1}, rate=k_creat1, restrict_to=model.EDGE1)
    R2 = spatialpy.Reaction(reactants=None, products={B: 1}, rate=k_creat2, restrict_to=model.EDGE2)
    R3 = spatialpy.Reaction(reactants={A: 1, B: 1}, products=None, rate=k_react)
    model.add_reaction([R1, R2, R3])

    # Define simulation timespan
    # self.set_timesteps(output_interval=1, num_steps=200, timestep_size=1)
    tspan = spatialpy.TimeSpan(range(500), timestep_size=1)
    model.timespan(tspan)
