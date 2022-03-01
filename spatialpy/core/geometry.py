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
from spatialpy.core.spatialpyerror import GeometryError

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
    only works for meshes that define triangles and tetrahedrons.
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
    Only works for meshes that define triangles and tetrahedrons.
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
