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

from spatialpy.core.spatialpyerror import LatticeError


class Lattice:
	"""
    Lattice class provides a method for creating parts of the spatial domain.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :raises LatticeError: if center is not and list, doesn't contain 3 values, or any value is not a float.
    """
    def __init__(self, center):
    	if isinstance(center, tuple):
    		center = list(center)
    	try:
    		self.center = [float(val) for val in center]
    		self.validate()
    	except ValueError as err:
    		raise LatticeError("Values in center must be of type float.") from err

    def apply(self, domain, geometry, transform=None, **kwargs):
    	"""
    	Fill a domain with particles within the lattice and geometry.

    	:param domain: Domain particles are to be added to.
    	:type domain: spatialpy.Domain

    	:param geometry: Geometry defining the region within the lattice for the particles
    	:type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry | spatialpy.Transformation

    	:param transform: Transformation function applied to each particle.
    	:type transform: function

    	:param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.
    	"""
    	raise LatticeError("Subclasses of spatialpy.Lattice must implement the apply() method")

    def validate(self):
    	"""
    	Validates the center point.
    	"""
    	if not isinstance(self.center, list):
    		raise LatticeError("center must be of type list.")
    	if len(self.center) < 3 or len(self.center) > 3:
    		raise LatticeError("center must be of length 3.")
