'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#!/usr/bin/env python3

import pickle
import unittest
import string
import spatialpy


#class diffusion_debug(spatialpy.Model):
#
#    def __init__(self, model_name="diffusion_debug_test"):
#        spatialpy.Model.__init__(self, model_name)
#
#        D_const = 0.01
#
#        A = spatialpy.Species(name="A", diffusion_constant=D_const)
#        self.add_species([A])
#
#        self.domain = spatialpy.Domain.create_2D_domain(
#            xlim=[-1, 1], ylim=[-1, 1], nx=50, ny=50, type_id=1.0,
#            mass=1.0, nu=1.0, fixed=True,  rho0=1.0, c0=1.0, P0=1.0
#        )
#
#        self.add_initial_condition(spatialpy.PlaceInitialCondition(A, 100000, [0,0,0]))
#
#        self.timestep_size=.1
#        self.num_timesteps=10
#        self.output_freq=1
#

# class testPeriodicDiffusion(spatialpy.Model):
#     def __init__(self, model_name="test1D"):
#         spatialpy.Model.__init__(self, model_name)
#         X = self.add_species(spatialpy.Species(name="X",  diffusion_constant=0.001))
#         self.domain = spatialpy.Domain.generate_unit_interval_mesh(nx=100, periodic=True)
#         self.add_initial_condition(spatialpy.PlaceInitialCondition(X, 1000))
#         #self.set_initial_condition_place_near({X:1000}, 0.1)
#         self.timespan(range(10))


class TestModelFunctionality(unittest.TestCase):

    def setUp(self):
        pass

    def test_single_letter_species_names(self):
        """ Test that we can have any single letter species name.

        NOTE:  this test must include 'model.run()' to ensure the compliation
        namespace is working correctly.
        """
        model = spatialpy.Model()
        for x in string.ascii_letters:
            with self.subTest(name=x):
                if x in spatialpy.Model.reserved_names:
                    with self.assertRaises(spatialpy.ModelError):
                        model.add_species(spatialpy.Species(x,0)) 
                else:
                    model.add_species(spatialpy.Species(x,0))
        model.set_timesteps(output_interval=1,num_steps=1,timestep_size=1)
        model.add_domain(spatialpy.Domain.create_2D_domain([0,1],[0,1],2,2))

        result = model.run(debug_level=0) #this will fail with Exception in the names checking is not correct



    def test_reaction_init(self):
        """ Test that we can instantate a Reaction is all the supported ways. """
        model = spatialpy.Model()
        s1 = spatialpy.Species('s1', D=0)
        s2 = spatialpy.Species('s2', D=0)
        model.add_species([s1, s2])
        p = spatialpy.Parameter('p',0,1)
        model.add_parameter(p)
        r1 = spatialpy.Reaction(reactants={s1:1},products={s2:1},rate=p)
        r2 = spatialpy.Reaction(reactants={'s1':1},products={'s2':1},rate='p')
        r3 = spatialpy.Reaction(reactants={'s1':1},products={'s2':1},rate=0.5)
        r4 = spatialpy.Reaction(reactants={'s1':1},products={'s2':1},rate=5)
        model.add_reaction([r1,r2,r3,r4])
        self.assertEqual(r1.propensity_function, r2.propensity_function)



