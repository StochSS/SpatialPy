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
from spatialpy.core.spatialpyerror import ParameterError

class Parameter():
    """
    Model of a rate paramter.
    A parameter can be given as a String expression (function) or directly as a scalar value.
    If given a String expression, it should be evaluable in the namespace of a parent Model.

    :param name: Name of the Parameter.
    :type name: str

    :param expression: Mathematical expression of Parameter.
    :type expression: str

    :param value: Parameter as value rather than expression.
    :type value: float
    """
    def __init__(self, name=None, expression=None):
        if name is None:
            raise ParameterError("name is required for a Parameter.")
        if not isinstance(name, str):
            raise ParameterError("Parameter name must be a string.")

        if expression is None:
            raise ParameterError("expression is required for a Parameter.")

        self.name = name
        self.value = None
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        if isinstance(expression, (int, float)):
            self.expression = str(expression)
        else:
            self.expression = expression

    def __str__(self):
        print_string = f"{self.name}: {str(self.expression)}"
        return print_string

    def _evaluate(self, namespace=None):
        """
        Evaluate the expression and return the (scalar) value.

        :param namespace: A dictionary containing key,value pairings of expressions and evaluable executions.
        :type namespace: dict
        """
        if namespace is None:
            namespace = {}
        try:
            self.value = (float(eval(self.expression, namespace)))
        except Exception as err:
            message = f"Could not evaluate expression '{self.expression}': {err}."
            raise ParameterError(message) from err
