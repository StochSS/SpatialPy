#!/usr/bin/env python

# import pickle
from spatialpy.Model import ModelError
import unittest
import spatialpy

# class SimpleDiffusion(spatialpy.Model):
#     """ Initial condition is a delta function at the center voxel.
#         The solution should be a Gaussian, up to the point where
#         the BC becomes important. """

#     def __init__(self, model_name="simple_diffusion"):
#         spatialpy.Model.__init__(self, model_name)
#         A = self.add_species(spatialpy.Species(
#             name="A", diffusion_constant=0.01))
#         # A unit square

#         # System constants
#         nxF, nyF = 50, 50     # number of particles in x and y-direction
#         L = 1         # characteristic lenght of the cavity (= width = height)
#         nW = 3         # number of wall points
#         rho0 = 1         # reference fluid density

#         # Compute domain bounds (including the boundary)
#         dx, dy = L/(nxF-1), L/(nyF-1)
#         xLim = ((0-(nW-1)*dx), 1+(nW-1)*dx)
#         yLim = ((0-(nW-1)*dy), 1+(nW-1)*dy)

#         # Discretization
#         # total number of particles in x-direction (including walls)
#         nxTot = nxF + 2*(nW-1)
#         # total number of particles in y-direction (including walls)
#         nyTot = nyF + 2*(nW-1)

#         # Compute volume and mass per particle
#         # in 2D simulations, consider z-lenght = 1
#         vol = (xLim[1]-xLim[0])*(yLim[1]-yLim[0])*1.0
#         # density * total volume / total number of particles
#         mPP = rho0*vol/(nxTot*nyTot)

#         self.mesh = spatialpy.Mesh.create_2D_domain(
#             xlim=xLim, ylim=yLim, nx=nxTot, ny=nyTot, type_id=1, mass=mPP, rho0=rho0, fixed=True)
#         # Place the A molecules in the voxel nearest the center of the square
#         # self.set_initial_condition_place_near({A:10000},point=[0.5,0.5])
#         self.add_initial_condition(spatialpy.ScatterInitialCondition(A, 10000))
#         self.timestep_size = 1e-4
#         self.timespan(numpy.linspace(0, 5, 200))


class diffusion_debug(spatialpy.Model):

    def __init__(self, model_name="diffusion_debug_test", diffusion_constant=0.01):
        spatialpy.Model.__init__(self, model_name)

        A = spatialpy.Species(name="A", diffusion_constant=diffusion_constant)
        self.add_species([A])

        self.mesh = spatialpy.Mesh.create_2D_domain(
            xlim=[-1, 1], ylim=[-1, 1], nx=50, ny=50, type_id=1.0,
            mass=1.0, nu=1.0, fixed=True,  rho0=1.0, c0=1.0, P0=1.0
        )

        self.add_initial_condition(spatialpy.PlaceInitialCondition(A, 1000, [0,0,0]))

        self.timestep_size=0.1
        self.num_timesteps=10
        self.output_freq=1

# class testPeriodicDiffusion(spatialpy.Model):
#     def __init__(self, model_name="test1D"):
#         spatialpy.Model.__init__(self, model_name)
#         X = self.add_species(spatialpy.Species(name="X",  diffusion_constant=0.001))
#         self.mesh = spatialpy.Mesh.generate_unit_interval_mesh(nx=100, periodic=True)
#         self.add_initial_condition(spatialpy.PlaceInitialCondition(X, 1000))
#         #self.set_initial_condition_place_near({X:1000}, 0.1)
#         self.timespan(range(10))


class TestSolverFunctionality(unittest.TestCase):

    def setUp(self):
        self.model = diffusion_debug()
        #self.periodic_model = testPeriodicDiffusion()

    def test_solver_io(self):
        """ Test that the initial value in the solver output file is the same as the input initial value. """
        model = diffusion_debug()
        result = model.run()
        A = result.get_species("A", 0)
        self.assertFalse((A-model.u0).any())

    # def test_zero_diffusion(self):
    #     """ Test that nothing happens if the diffusion is set to zero. """
    #     model = diffusion_debug(diffusion_constant=0.0)
    #     result = model.run()
    #     A = result.get_species("A", -1)
    #     self.assertFalse((A - model.u0).any())

    # def test_same_seed(self):
    #     """ Test that the output is the same if the same seed is used, edxplicit solver creation  """
    #     solver = spatialpy.Solver(self.model)
    #     result1 = solver.run(seed=1)
    #     result2 = solver.run(seed=1)
    #     self.assertEqual(result1, result2)

    # def test_same_seed2(self):
    #     """ Test that the output is the same if the same seed is used, use model.run()  """
    #     result1 = self.model.run(seed=1)
    #     result2 = self.model.run(seed=1)
    #     self.assertEqual(result1, result2)

    # def test_different_seeds(self):
    #     """ Test that the output is different if different seeds are used. """
    #     solver = spatialpy.Solver(self.model)
    #     result1 = solver.run(seed=1)
    #     result2 = solver.run(seed=100)
    #     self.assertNotEqual(result1, result2)

    # def test_default_seed(self):
    #     """ Test that the output is different if no seed is given (default set on C level). """
    #     solver = spatialpy.Solver(self.model)
    #     result1 = solver.run()
    #     result2 = solver.run()
    #     self.assertNotEqual(result1, result2)

    # def test_mesh_pickle(self):
    #     meshstr = pickle.dumps(self.model.mesh)
    #     mesh = pickle.loads(meshstr)

    # def test_model_pickle(self):
    #     """ Test that the model is picklable. We do not compare models directly, but rather the results after simulation. """
    #     model = self.model
    #     model_str = pickle.dumps(model)
    #     model2 = pickle.loads(model_str)
    #     result1 = model.run(seed=1)
    #     result2 = model2.run(seed=1)
    #     self.assertEqual(result1,result2)

    # def test_solver_pickle(self):
    #     """ Test that the model, solver and result objects are pickleable. """
    #     sol = pyurdme.nsmsolver.NSMSolver(self.model)
    #     sol_str = pickle.dumps(sol)
    #     sol2 = pickle.loads(sol_str)
    #     result1 = sol.run(seed=1)
    #     result2 = sol2.run(seed=1)
    #     self.assertEqual(result1,result2)

    # def test_result_pickle(self):
    #     """ Test that the result object is picklable. """
    #     sol = pyurdme.nsmsolver.NSMSolver(self.model)
    #     result = sol.run(seed=1)
    #     result_str = pickle.dumps(result)
    #     result2 = pickle.loads(result_str)
    #     self.assertEqual(result,result2)

    def test_run_ensemble(self):
        """ Test the running of ensembles of runs """
        result_list = self.model.run(3)
        self.assertEqual(len(result_list), 3)

    # def test_1D_periodic_boundary(self):
    #     """ Test if periodic boundary conditions are working. """
    #     result = self.periodic_model.run()
    #     self.assertTrue(result.get_species("X", timepoints=1)[-1] > 0)

    # def test_1D_periodic_boundary_pickle(self):
    #     """ Test if periodic boundary conditions are working. """
    #     model_str = pickle.dumps(self.periodic_model)
    #     model2 = pickle.loads(model_str)
    #     self.periodic_model.assemble()
    #     model2.assemble()
    #     self.assertEqual(self.periodic_model.mesh.get_num_dof_voxels(), model2.mesh.get_num_dof_voxels())


if __name__ == '__main__':
    unittest.main()
