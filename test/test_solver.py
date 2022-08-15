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

#!/usr/bin/env python3

import os
import tempfile
import subprocess
import pickle
import unittest
import numpy
import spatialpy
from spatialpy.solvers.build_expression import BuildExpression, ExpressionConverter

def create_diffusion_debug(model_name="diffusion_debug_test", parameter_values=None):
    model = spatialpy.Model(model_name)

    D_const = 0.01

    A = spatialpy.Species(name="A", diffusion_coefficient=D_const)
    model.add_species([A])

    domain = spatialpy.Domain.create_2D_domain(
        xlim=[-1, 1], ylim=[-1, 1], numx=50, numy=50, type_id=1,
        mass=1.0, nu=1.0, fixed=True,  rho0=1.0, c0=1.0, P0=1.0
    )
    model.add_domain(domain)

    model.add_initial_condition(spatialpy.PlaceInitialCondition(A, 100000, [0,0,0]))

    tspan = spatialpy.TimeSpan.linspace(t=10, num_points=11, timestep_size=0.1)
    model.timespan(tspan)
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

class ExpressionTestCase:
    """
    Each test expression consists of a dict of argument names, an expression, and a list of
    values to be passed as arguments to the given expression.
    """
    def __init__(self, args: "dict[str, str]", expression: "str", values: "list[list[float]]"):
        self.args = args
        self.expression = expression
        self.values = values

class TestSolverFunctionality(unittest.TestCase):

    expressions = [
        # Each test expression consists of a dict of args, an expression string, and a list of arg values.
        # Asserts that single operations work.
        ExpressionTestCase({"x": "x"}, "x*2", [
            [0.0], [1.0], [-1.0], [9.999], [-9.999],
        ]),
        # Asserts that order of operations is being evaluated properly.
        ExpressionTestCase({"x": "x"}, "x*2 + x/2 - (x*3)^2 + x/3^2", [
            [0.0], [1.0], [-1.0], [3.333], [-3.333], [9.8765], [-9.8765]
        ]),
        # Asserts that order of operations is evaluated properly with multiple variables.
        ExpressionTestCase({"x": "x", "y": "y"}, "(x-1)*y^2+x", [
            [1.0, 2.4], [5.1, 0.0], [5.1, 1.0], [5.1, -1.0], [9.8765, -1.0], [-1.0, 9.8765],
        ]),
        # Asserts complex order of operations with a large number of variables.
        ExpressionTestCase({"x": "x", "y": "y", "z": "z"}, "(x^2/y^2/z^2)/x^2/y^2/z^2**1/x**1/y**1/z", [
            [5.1, 0.1, 2.0], [0.1, 5.1, 2.0], [2.0, 0.1, 5.1], [2.0, 5.1, 0.1],
        ]),
    ]
    comparisons = [
        # Asserts that single comparison expressions work.
        ExpressionTestCase({"x": "x"}, "x > 0", [
            [100], [0], [0.001], [-1],
        ]),
        ExpressionTestCase({"x": "x", "y": "y"}, "x > y", [
            [100, 99], [99, 100], [-10, 10], [10, -10],
            [0.001, 0.0], [0.0, 0.001], [-99.999, -99.998]
        ]),
        # Asserts that single boolean operators work.
        ExpressionTestCase({"x": "x", "y": "y"}, "x > 0 and y < x", [
            [100, 99], [99, 100], [0, -100], [-0.001, -99.0], [0, 0.001], [-0.001, 0]
        ]),
        # Asserts that nested boolean operators work.
        ExpressionTestCase({"x": "x", "y": "y"}, "x > 0 and y < 10 and x > y", [
            [100, 9], [0.01, 0.00], [100, 200], [0.01, 0.02],
            [0, 0], [-0.01, -0.02], [-0.01, 0],
        ]),
        # Asserts that both && and || work.
        ExpressionTestCase({"x": "x", "y": "y"}, "x > 0 and y < 10 or y > 100", [
            [10, 9], [0.01, 9.99], [0, 10], [-1.0, -1.0],
        ]),
        # Asserts that nested boolean operators properly respect order of operations.
        ExpressionTestCase({"x": "x", "y": "y", "z": "z"}, "x^2>x and y<y^2 or z^2!=z^3 and y!=z", [
            [1.0, 1.0, 1.0], [99.9, 99.9, 100.0],
            [0.0, -1.0, 99.9], [-1.0, -1.0, 0.00],
        ]),
    ]

    def setUp(self):
        self.model = create_diffusion_debug()
        # self.periodic_model = create_periodic_diffusion()

    def test_solver_io(self):
        """ Test that the initial value in the solver output file is the same as the input initial value. """
        result = self.model.run()
        A = result.get_species("A", 0)
        self.assertFalse((A - self.model.u0).any())

    # def test_zero_diffusion(self):
    #     """ Test that nothing happens if the diffusion is set to zero. """
    #     self.model.listOfSpecies['A'].diffusion_coefficient = 0.0
    #     result = self.model.run()
    #     A = result.get_species("A", -1)
    #     self.assertFalse((A - self.model.u0).any())

    def test_same_seed(self):
        """ Test that the output is the same if the same seed is used, explicit solver creation  """
        solver = spatialpy.Solver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_same_seed2(self):
        """ Test that the output is the same if the same seed is used, use model.run()  """
        result1 = self.model.run(seed=1)
        result2 = self.model.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_different_seeds(self):
        """ Test that the output is different if different seeds are used. """
        solver = spatialpy.Solver(self.model)
        result1 = solver.run(seed=1)
        result2 = solver.run(seed=100)
        self.assertFalse(result1 == result2)

    def test_default_seed(self):
        """ Test that the output is different if no seed is given (default set on C level). """
        solver = spatialpy.Solver(self.model)
        result1 = solver.run()
        result2 = solver.run()
        self.assertFalse(result1 == result2)

    def test_mesh_pickle(self):
        meshstr = pickle.dumps(self.model.domain)
        domain = pickle.loads(meshstr)

    def test_model_pickle(self):
        """ Test that the model is picklable. We do not compare models directly, but rather the results after simulation. """
        model = self.model
        model_str = pickle.dumps(model)
        model2 = pickle.loads(model_str)
        result1 = model.run(seed=1)
        result2 = model2.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_solver_pickle(self):
        """ Test that the model, solver and result objects are pickleable. """
        sol = spatialpy.Solver(self.model)
        sol_str = pickle.dumps(sol)
        sol2 = pickle.loads(sol_str)
        result1 = sol.run(seed=1)
        result2 = sol2.run(seed=1)
        self.assertTrue(result1 == result2)

    def test_result_pickle(self):
        """ Test that the result object is picklable. """
        sol = spatialpy.Solver(self.model)
        result = sol.run(seed=1)
        result_str = pickle.dumps(result)
        result2 = pickle.loads(result_str)
        self.assertTrue(result == result2)

    def test_run_ensemble(self):
        """ Test the running of ensembles of runs """
        result_list = self.model.run(3)
        self.assertEqual(len(result_list), 3)

    def test_solver_expressions(self):
        """
        Ensure that expression conversions to C++ result in (roughly) equivalent values as Python.
        """
        tmpdir = tempfile.mkdtemp()
        src_path = os.path.join(os.path.dirname(__file__), "assets", "evaluate.c")
        exe_path = os.path.join(tmpdir, "test")

        def build(expr_args: "list[str]", expr_str: "str", use_bool=False):
            args = ["gcc", "-o", exe_path, src_path, "-lm"]
            expr_num = str(len(expr_args))
            expr_args = ",".join(expr_args)
            args.append(f"-DEXP{expr_num}({expr_args})=({expr_str})")
            if use_bool:
                args.append("-DUSE_BOOLEAN")
            subprocess.check_call(args)

        def run(args: "list[str]") -> str:
            args.insert(0, exe_path)
            stdout = subprocess.check_output(args)
            return stdout.decode("ascii")

        def test_expressions(expressions: "list[ExpressionTestCase]", use_bool=False):
            for entry in expressions:
                expression = ExpressionConverter.convert_str(entry.expression)
                expr = BuildExpression(namespace=entry.args)
                cpp_expr = expr.getexpr_cpp(expression)
                with self.subTest(msg="Evaluating converted C expressions",
                                  expression=entry.expression,
                                  c_expression=cpp_expr):
                    py_args = ",".join(entry.args.keys())
                    py_func = eval(f"lambda {py_args}: {expression}")

                    for value_set in entry.values:
                        value_str = [str(val) for val in value_set]
                        with self.subTest(values=",".join(value_str)):
                            expect = py_func(*value_set)
                            build(list(entry.args.values()), cpp_expr, use_bool)
                            if use_bool:
                                result_cpp = bool(int(run(value_str)))
                                self.assertTrue(expect == result_cpp)
                            else:
                                result_cpp = float(run(value_str))
                                self.assertAlmostEqual(expect, result_cpp, places=3)

        try:
            # Test expressions which return a float value
            test_expressions(self.expressions, False)
            # Test boolean and comparator expressions
            test_expressions(self.comparisons, True)

        finally:
            import shutil
            shutil.rmtree(tmpdir)


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
    #     self.assertEqual(self.periodic_model.domain.get_num_dof_voxels(), model2.domain.get_num_dof_voxels())


if __name__ == '__main__':
    unittest.main()
