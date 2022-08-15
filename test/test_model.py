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

import pickle
import unittest
import string
import spatialpy

def create_diffusion_debug(model_name="diffusion_debug_test", parameter_values=None):
    model = spatialpy.Model(model_name)

    D_const = 0.01

    A = spatialpy.Species(name="A", diffusion_coefficient=D_const)
    model.add_species([A])

    domain = spatialpy.Domain.create_2D_domain(
       xlim=[-1, 1], ylim=[-1, 1], nx=50, ny=50, type_id=1.0,
       mass=1.0, nu=1.0, fixed=True,  rho0=1.0, c0=1.0, P0=1.0
    )
    model.add_domain(domain)

    model.add_initial_condition(spatialpy.PlaceInitialCondition(A, 100000, [0, 0, 0]))

    model.set_timesteps(output_interval=1, num_steps=10, timestep_size=0.1)
    return model

def create_periodic_diffusion(model_name="test1D", parameter_values=None):
    model = spatialpy.Model(model_name)
    
    X = spatialpy.Species(name="X",  diffusion_coefficient=0.001)
    model.add_species(X)

    domain = spatialpy.Domain.generate_unit_interval_mesh(nx=100, periodic=True)
    model.add_domain(domain)

    model.add_initial_condition(spatialpy.PlaceInitialCondition(X, 1000))
    
    tspan = spatialpy.TimeSpan(range(10))
    model.timespan(tspan)
    return model

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


    def test_data_function(self):
        """ Test if the data function is working correctly. """
        model = spatialpy.Model()
        for x in spatialpy.Model.reserved_names:
            with self.subTest(name=x):
                with self.assertRaises(spatialpy.ModelError):
                    df = spatialpy.DataFunction(name=x)
                    model.add_data_function(df)
        # check if it actually changes at different places in the domain
        df = spatialpy.DataFunction(name="df")
        df.map = lambda x: x[0]*10000
        model.add_data_function(df)
        model.set_timesteps(output_interval=1,num_steps=1,timestep_size=1)
        model.add_domain(spatialpy.Domain.create_2D_domain([0,1],[0,1],2,2))
        model.add_species(spatialpy.Species('A',0))
        model.add_reaction(spatialpy.Reaction(products={'A':1},propensity_function="df"))
        df_res = model.run()
        traj = df_res.get_species('A',-1)
        self.assertTrue(traj[0] == 0)
        self.assertTrue(traj[1] == 0)
        self.assertTrue(traj[2] > 0)
        self.assertTrue(traj[3] > 0)




    def test_reaction_init(self):
        """ Test that we can instantate a Reaction is all the supported ways. """
        model = spatialpy.Model()
        s1 = spatialpy.Species('s1', diffusion_coefficient=0)
        s2 = spatialpy.Species('s2', diffusion_coefficient=0)
        model.add_species([s1, s2])
        p = spatialpy.Parameter('p',1)
        model.add_parameter(p)
        r1 = spatialpy.Reaction(reactants={s1:1},products={s2:1},rate=p)
        r2 = spatialpy.Reaction(reactants={'s1':1},products={'s2':1},rate='p')
        r3 = spatialpy.Reaction(reactants={'s1':1},products={'s2':1},rate=0.5)
        r4 = spatialpy.Reaction(reactants={'s1':1},products={'s2':1},rate=5)
        model.add_reaction([r1,r2,r3,r4])
        self.assertEqual(r1.propensity_function, r2.propensity_function)

    def test_parameters(self):
        """ Test that we can create and add Parameters to a model."""
        m1 = spatialpy.Model()
        m1.set_timesteps(output_interval=1,num_steps=1,timestep_size=1)
        m1.add_domain(spatialpy.Domain.create_2D_domain([0,1],[0,1],2,2))

        p1 = spatialpy.Parameter(name="p1",expression="1+1")
        m1.add_parameter(p1)
        p2 = spatialpy.Parameter(name="p2",expression="2")
        p3 = spatialpy.Parameter(name="p3",expression="2")
        m1.add_parameter([p2,p3])

        s1 = spatialpy.Solver(m1)
        s1.compile()













