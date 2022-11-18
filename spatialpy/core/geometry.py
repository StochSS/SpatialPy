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

from spatialpy.core.spatialpyerror import GeometryError

class CombinatoryGeometry:
    """
    Combinatory Geometry class uses one or more Geometry class for inclusion or exclusion of
    multiple geometric shapes.

    :param formula: Boolean logic formula to be evaluated. Example use: "(geo1 or geo2) and not geo3".
    :type formula: str

    :param geo_namespace: Namespace used when evaluating the formula. Example use:
        {'geo1': Geometry1(), 'geo2': Geometry2(), 'geo3': Geometry3()} where 'geo1',
        'geo2', and 'geo3' are found in the formula.
    :type geo_namespace: dict
    """
    def __init__(self, formula, geo_namespace):
        self.formula = formula
        self.geo_namespace = geo_namespace

        self.validate()

    def inside(self, point, on_boundary):
        """
        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: True if the particle satisfies boolean condition for all geometies, else False.
        :rtype: bool

        :raises GeometryError: If any geometries inside method does not return a bool.
        """
        namespace = {}
        for name, geometry in self.geo_namespace.items():
            val = geometry.inside(point, on_boundary)
            if not isinstance(val, (bool, numpy.bool_)):
                errmsg = f"{name} is not a valid Geometry obj. Reason given: inside() method must return a bool"
                raise GeometryError(errmsg)
            namespace[name] = val
        return eval(self.formula, {}, namespace)

    def validate(self):
        """
        Validate the combinatory geometry object.
        """
        keys = list(self.geo_namespace.keys())
        for name in keys:
            if name not in self.formula:
                raise GeometryError("geo_namespace entries must be in the formula.")
        keys.extend(['and', 'or', 'not', '(', ')'])
        formula = self.formula
        for item in keys:
            formula = formula.replace(item, "")
        if formula.strip() != "":
            errmsg = f"formula may only contain boolean operators or geometries in geo_namespace not {formula.split(' ')}"
            raise GeometryError(errmsg)

class Geometry:
    """
    Geometry class provides a method for tagging parts of the spatial domain as separate parts.
    """
    def __init__(self):
        pass

    def inside(self, point, on_boundary):
        """
        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: True if the particle is in the geometric shape, else False.
        :rtype: bool
        """
        raise GeometryError("Subclasses of spatialpy.Geometry must implement the inside() method")

class GeometryAll(Geometry):
    """
    Mark all particles.
    """
    def inside(self, point, on_boundary):
        """
        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: True if the particle is in the domain.
        :rtype: bool
        """
        return True

class GeometryExterior(Geometry):
    """
    Mark particles that are on the edge of the domain.
    only works for domains that define triangles and tetrahedrons.
    """
    def inside(self, point, on_boundary):
        """
        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: True if the particle is on the edge of the domain, else False.
        :rtype: bool
        """
        return on_boundary

class GeometryInterior(Geometry):
    """
    Mark particles that are not on the edge of the domain.
    Only works for domains that define triangles and tetrahedrons.
    """
    def inside(self, point, on_boundary):
        """
        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: True if the particle is not on the edge of the domain, else False.
        :rtype: bool
        """
        return not on_boundary
