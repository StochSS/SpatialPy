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

from spatialpy.core.domain import Domain
from spatialpy.core.geometry import Geometry, CombinatoryGeometry
from spatialpy.core.lattice import Lattice
from spatialpy.core.spatialpyerror import TransformationError


class Transformation:
    """
    Transformation class provides a methods for applying and reversing
        transformations to parts the spatial domain.

    :param geometry: Geometry classed used once the transformation is reversed.
        Transformation.inside wraps Geometry.inside.
    :type geometry: spatialpy.Geometry

    :param lattice: Lattice classed used when applying transformations.
        Transformation.apply wraps Lattice.apply.
    :type lattice: spatialpy.Lattice

    :param transformation: Transformation to be applied after applying this
        transformation.
    :type transformation: spatialpy.Transformation

    :raises TransformationError: if the provided Geometry, Lattice, or
        Transformation are invalid.
    """
    def __init__(self, geometry=None, lattice=None, transformation=None, skip_validate=False):
        self.geometry = geometry
        self.lattice = lattice
        self.transformation = transformation

        if not skip_validate:
            self.validate()

    def apply(self, domain, geometry, **kwargs):
        """
        Wrapper for Lattice.apply that passes in a transformation function used to transform
            particles within the region defined by the lattice and geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param geometry: Geometry defining the region within the lattice in
            which particles are restricted to.
        :type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.

        :returns: The results of Lattice.apply.
        :rtype: int

        :raises TransformationError: if a lattice was not provided.
        """
        if self.lattice is None:
            raise TransformationError("Transformations must wrap a lattice to use the apply function.")
            
        return self.lattice.apply(domain, geometry, transform=self.transform, **kwargs)

    def inside(self, point, on_boundary):
        """
        Wrapper for Geometry.inside.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: The results of Geometry.inside.
        :rtype: bool

        :raises TransformationError: if a geometry was not provided.
        """
        if self.geometry is None:
            raise TransformationError("Transformations must wrap a geometry to use the inside function.")
            
        r_point = self.reverse_transform(point)
        return self.geometry.inside(r_point, on_boundary)

    def reverse_transform(self, point):
        """
        Reverses the define transformation for the given point.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :returns: The point prior to any transformations.
        :rtype: float[3]
        """
        raise TransformationError(
            "Subclasses of spatialpy.Transformation must implement the reverse_transform() method."
        )

    def transform(self, point):
        """
        Applies the define transformation to the given point.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :returns: The point prior to any transformations.
        :rtype: float[3]
        """
        raise TransformationError(
            "Subclasses of spatialpy.Transformation must implement the transform() method."
        )

    def validate(self):
        """
        Validate the geometry, lattice, and transformation attributes.
        """
        if self.geometry is not None and \
            not (isinstance(self.geometry, Geometry) or type(self.geometry).__name__ == "Geometry"):
            raise TransformationError("geometry must be of type spatialpy.Geometry.")

        if self.lattice is not None and \
            not (isinstance(self.lattice, Lattice) or type(self.lattice).__name__ == "Lattice"):
            raise TransformationError("lattice must be of type spatialpy.Lattice.")

        if self.transformation is not None and \
            not (isinstance(self.transformation, Transformation) or type(self.transformation).__name__ == "Transformation"):
            raise TransformationError("transformation must be of type spatialpy.Transformation.")

class TranslationTransformation(Transformation):
    """
    Translation transformation class provides a methods for applying and reversing
        translation transformations by an arbitrary vector to parts the spatial
        domain.

    :param vector: Arbitrary vector by which the region is translated by.
        Basic example: [[2, 2, 3], [4, 5, 5]]
    :type vector: list

    :param geometry: Geometry classed used once the transformation is reversed.
        Transformation.inside wraps Geometry.inside.
    :type geometry: spatialpy.Geometry

    :param lattice: Lattice classed used when applying transformations.
        Transformation.apply wraps Lattice.apply.
    :type lattice: spatialpy.Lattice

    :param transformation: Transformation to be applied after applying this
        transformation.
    :type transformation: spatialpy.Transformation

    :raises TransformationError: if the provided Geometry, Lattice, or
        Transformation are invalid.
    """
    def __init__(self, vector, geometry=None, lattice=None, transformation=None):
        super().__init__(
            geometry=geometry, lattice=lattice, transformation=transformation, skip_validate=True
        )

        if isinstance(vector, (list, tuple)):
            vector = numpy.array(vector)
        elif not isinstance(vector, numpy.ndarray):
            raise TransformationError("vector must be of type list.")
        
        try:
            if vector.shape == (3,):
                o_vector = vector
            elif vector.shape == (1, 3):
                o_vector = vector[0]
            elif vector.shape == (2, 3):
                if numpy.count_nonzero(vector[0]) == 0:
                    o_vector = vector[1]
                else:
                    o_vector = vector[1] - vector[0]
        except TypeError:
            raise TransformationError("vector coordinates must be of type float.")
        
        self.vector = o_vector

        self.validate()

    # def __execute(self, point, vector):
    
    def reverse_transform(self, point):
        """
        Reverses the defined rotation transformation for the given point.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :returns: The point prior to any transformations.
        :rtype: float[3]
        """
        if self.transformation is not None:
            point = self.transformation.reverse_transform(point)
        
        sin_t = numpy.sin(-self.theta)
        cos_t = numpy.cos(-self.theta)
        
        return self.__execute(point, sin_t, cos_t)

    def transform(self, point):
        """
        Applies the defined rotation transformation to the given point.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :returns: The point prior to any transformations.
        :rtype: float[3]
        """
        sin_t = numpy.sin(self.theta)
        cos_t = numpy.cos(self.theta)
        
        t_point = self.__execute(point, sin_t, cos_t)
        
        if self.transformation is not None:
            return self.transformation.transform(t_point)
        
        return t_point

    def validate(self):
        """
        Validate the rotation transformation attributes.
        """
        super().validate()

        if not isinstance(self.theta, (int, float)):
            raise TransformationError("angle must be of type float.")

class RotationTransformation(Transformation):
    """
    Rotation transformation class provides a methods for applying and reversing
        rotation transformations around an arbitrary vector by a given angle to
        parts the spatial domain.

    :param vector: Arbitrary vector by which the region is rotated around.
        Basic example: [[2, 2, 3], [4, 5, 5]]
    :type vector: list

    :param angle: Angle by which to rotate the region.
    :type angle: float

    :param geometry: Geometry classed used once the transformation is reversed.
        Transformation.inside wraps Geometry.inside.
    :type geometry: spatialpy.Geometry

    :param lattice: Lattice classed used when applying transformations.
        Transformation.apply wraps Lattice.apply.
    :type lattice: spatialpy.Lattice

    :param transformation: Transformation to be applied after applying this
        transformation.
    :type transformation: spatialpy.Transformation

    :raises TransformationError: if the provided Geometry, Lattice, or
        Transformation are invalid.
    """
    def __init__(self, vector, angle, geometry=None, lattice=None, transformation=None):
        super().__init__(
            geometry=geometry, lattice=lattice, transformation=transformation, skip_validate=True
        )

        if isinstance(vector, (list, tuple)):
            vector = numpy.array(vector)
        elif not isinstance(vector, numpy.ndarray):
            raise TransformationError("vector must be of type list.")
        
        try:
            if vector.shape == (3,):
                o_vector = vector
            elif vector.shape == (1, 3):
                o_vector = vector[0]
            elif vector.shape == (2, 3):
                if numpy.count_nonzero(vector[0]) == 0:
                    o_vector = vector[1]
                else:
                    o_vector = vector[1] - vector[0]
                    
            magnitude = numpy.sqrt(sum(o_vector**2))
            if magnitude == 1.0:
                u_vector = o_vector
            else:
                u_vector = o_vector / magnitude
        except TypeError:
            raise TransformationError("vector coordinates must be of type float.")
        
        self.theta = angle
        self.vector = u_vector

        self.validate()

    def __execute(self, point, sin_t, cos_t):
        x_y = self.vector[0] * self.vector[1] * (1 - cos_t)
        x_z = self.vector[0] * self.vector[2] * (1 - cos_t)
        y_z = self.vector[1] * self.vector[2] * (1 - cos_t)
        x_sin_t = self.vector[0] * sin_t
        y_sin_t = self.vector[1] * sin_t
        z_sin_t = self.vector[2] * sin_t
        
        matrix = [
            [cos_t + self.vector[0]**2 * (1 - cos_t), x_y - z_sin_t, x_z + y_sin_t],
            [x_y + z_sin_t, cos_t + self.vector[1]**2 * (1 - cos_t), y_z - x_sin_t],
            [x_z - y_sin_t, y_z + x_sin_t, cos_t + self.vector[2]**2 * (1 - cos_t)]
        ]
        
        return numpy.matmul(matrix, point)
    
    def reverse_transform(self, point):
        """
        Reverses the defined rotation transformation for the given point.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :returns: The point prior to any transformations.
        :rtype: float[3]
        """
        if self.transformation is not None:
            point = self.transformation.reverse_transform(point)
        
        sin_t = numpy.sin(-self.theta)
        cos_t = numpy.cos(-self.theta)
        
        return self.__execute(point, sin_t, cos_t)

    def transform(self, point):
        """
        Applies the defined rotation transformation to the given point.

        :param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

        :returns: The point prior to any transformations.
        :rtype: float[3]
        """
        sin_t = numpy.sin(self.theta)
        cos_t = numpy.cos(self.theta)
        
        t_point = self.__execute(point, sin_t, cos_t)
        
        if self.transformation is not None:
            return self.transformation.transform(t_point)
        
        return t_point

    def validate(self):
        """
        Validate the rotation transformation attributes.
        """
        super().validate()

        if not isinstance(self.theta, (int, float)):
            raise TransformationError("angle must be of type float.")
