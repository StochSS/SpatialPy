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
    	raise TransformationError(
    		"Subclasses of spatialpy.Transformation must implement the apply() method."
    	)

    def inside(self, point, on_boundary):
    	"""
    	Wrapper for Lattice.apply that passes in a transformation function used to transform
    		particles within the region defined by the lattice and geometry.

    	:param point: X, Y, Z coodinates for the particle.
        :type point: float[3]

    	:param on_boundary: Indicates whether a particle is on the edge of the domain.
        :type on_boundary: bool

        :returns: The results of Geometry.inside.
        :rtype: bool

        :raises TransformationError: if a geometry was not provided.
    	"""
    	raise TransformationError(
    		"Subclasses of spatialpy.Transformation must implement the inside() method."
    	)

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
