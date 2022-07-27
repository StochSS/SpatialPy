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

    :raises LatticeError: if center is not and list, doesn't contain 3 values,
    	or any value is not a float.
    """
    def __init__(self, center=None):
    	if center is None:
    		center = [0] * 3
    	elif isinstance(center, tuple):
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

    	:param geometry: Geometry defining the region within the lattice in
    		which particles are restricted to.
    	:type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry |
    		spatialpy.Transformation

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

class CartesianLattice(Lattice):
	"""
    Cartesian lattice class provides a method for creating parts of the spatial
    	domain within a cartesian coordinate system.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :param xmin: Minimum x value of the lattice.
    :type xmin: float

    :param xmax: Maximum x value of the lattice.
    :type xmax: float

    :param ymin: Minimum y value of the lattice (optional, defaults to xmin).
    :type ymin: float

    :param ymax: Maximum y value of the lattice (optional, defaults to xmax).
    :type ymax: float

    :param zmin: Minimum z value of the lattice (optional, defaults to xmin).
    :type zmin: float
    
    :param zmax: Maximum z value of the lattice (optional, defaults to xmax).
    :type zmax: float

    :param deltax: Distance between two particles on the X-Axis.
    :type deltax: float

    :param deltay: Distance between two particles on the Y-Axis (optional, defaults to deltax).
    :type deltay: float

    :param deltaz: Distance between two particles on the Z-Axis (optional, defaults to deltax).
    :type deltaz: float

    :raises LatticeError: if center is not and list, doesn't contain 3 values,
    	or any value is not a float.
    """
    def __init__(self, xmin, xmax, deltax, center=None, ymin=None,
    			 ymax=None, zmin=None, zmax=None, deltay=None, deltaz=None):
    	super().__init__(center)

    	if ymin is None:
    		ymin = xmin
    	if ymax is None:
    		ymax = xmax
    	if zmin is None:
    		zmin = xmin
    	if zmax is None:
    		zmax = xmax

    	self.xmin = xmin
    	self.xmax = xmax
    	self.ymin = ymin
    	self.ymax = ymax
    	self.zmin = zmin
    	self.zmax = zmax

    	if deltay is None:
    		deltay = deltax
    	if deltaz is None:
    		deltaz = deltax

    	self.deltax = deltax
    	self.deltay = deltay
    	self.deltaz = deltaz

    	self.validate()

    def apply():

    def validate():
    	"""
    	Validate the cartesian lattice dependencies.
    	"""
    	if not isinstance(self.xmin, (int, float)):
    		raise LatticeError("xmin must be of type float.")
		if not isinstance(self.xmax, (int, float)):
			raise LatticeError("xmax must be of type float.")
		if not isinstance(self.ymin, (int, float)):
			raise LatticeError("ymin must be of type float.")
		if not isinstance(self.ymax, (int, float)):
			raise LatticeError("ymax must be of type float.")
		if not isinstance(self.zmin, (int, float)):
			raise LatticeError("zmin must be of type float.")
		if not isinstance(self.zmax, (int, float)):
			raise LatticeError("zmax must be of type float.")
		if self.xmax <= self.xmin:
			raise LatticeError("xmax must be greater than xmin.")
		if self.ymax <= self.ymin:
			raise LatticeError("ymax must be greater than ymin.")
		if self.zmax <= self.zmin:
			raise LatticeError("zmax must be greater than zmin.")
		if not isinstance(deltax, (int, float)):
			raise LatticeError("deltax must be of type float.")
		if deltax <= 0:
			raise LatticeError("deltax must be greater than 0.")
		if not isinstance(deltay, (int, float)):
			raise LatticeError("deltay must be of type float.")
		if deltay <= 0:
			raise LatticeError("deltay must be greater than 0.")
		if not isinstance(deltaz, (int, float)):
			raise LatticeError("deltaz must be of type float.")
		if deltaz <= 0:
			raise LatticeError("deltaz must be greater than 0.")
