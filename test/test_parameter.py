'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2019 - 2022 SpatialPy developers.

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
import unittest

import spatialpy
from spatialpy import Parameter
from spatialpy import ParameterError

class TestParameter(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for spatialpy.Parameter.
    ################################################################################################
    '''
    def test_constructor(self):
        """ Test the Parameter constructor. """
        parameter = Parameter(name="test_parameter", expression="0.5")
        self.assertEqual(parameter.name, "test_parameter")
        self.assertEqual(parameter.expression, "0.5")


    def test_constructor__no_name(self):
        """ Test the Parameter constructor without name. """
        with self.assertRaises(ParameterError):
            parameter = Parameter(expression="0.5")


    def test_constructor__name_not_str(self):
        """ Test the Parameter constructor with non-str name. """
        with self.assertRaises(ParameterError):
            parameter = Parameter(name=0, expression="0.5")


    def test_constructor__no_expression(self):
        """ Test the Parameter constructor without expression. """
        with self.assertRaises(ParameterError):
            parameter = Parameter(name="test_parameter")


    def test_constructor__int_expression(self):
        """ Test the Parameter constructor with int expression. """
        parameter = Parameter(name="test_parameter", expression=1)
        self.assertEqual(parameter.expression, "1")


    def test_constructor__float_expression(self):
        """ Test the Parameter constructor with float expression. """
        parameter = Parameter(name="test_parameter", expression=0.5)
        self.assertEqual(parameter.expression, "0.5")


    def test___str__(self):
        """ Test Parameter.__str__ method. """
        parameter = Parameter(name="test_parameter", expression="0.5")
        self.assertIsInstance(str(parameter), str)


    def test__evaluate(self):
        """ Test Parameter._evaluate method. """
        parameter = Parameter(name="test_parameter", expression="0.5")
        parameter._evaluate()
        self.assertEqual(parameter.value, 0.5)


    def test__evaluate__parameter_in_namespace(self):
        """ Test Parameter._evaluate method with parameter in namespace. """
        parameter = Parameter(name="test_parameter", expression="k1 + 0.5")
        parameter._evaluate(namespace={"k1": 3})
        self.assertEqual(parameter.value, 3.5)


    def test__evaluate__species_in_namespace(self):
        """ Test Parameter._evaluate method with species in namespace. """
        parameter = Parameter(name="test_parameter", expression="S0 + 0.5")
        parameter._evaluate(namespace={"S0": 100})
        self.assertEqual(parameter.value, 100.5)


    def test__evaluate__improper_expression(self):
        """ Test Parameter._evaluate method with invalid expression. """
        parameter = Parameter(name="test_parameter", expression="[0.5]")
        with self.assertRaises(ParameterError):
            parameter._evaluate()


    def test__evaluate__param_not_in_namespace(self):
        """ Test Parameter._evaluate method with arg missing from namespace. """
        parameter = Parameter(name="test_parameter", expression="k1 + 0.5")
        with self.assertRaises(ParameterError):
            parameter._evaluate()
