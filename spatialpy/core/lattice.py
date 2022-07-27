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
from spatialpy.core.spatialpyerror import LatticeError


class Lattice:
    """
    Lattice class provides a method for creating parts of the spatial domain.

    :param center: The center point of the lattice (optional, defaults to [0, 0, 0]).
    :type center: float[3] | float(3)

    :raises LatticeError: if center is not and list, doesn't contain 3 values,
        or any value is not a float.
    """
    def __init__(self, center=None, skip_validate=False):
        if center is None:
            center = [0] * 3
        elif isinstance(center, tuple):
            center = list(center)
        try:
            self.center = [float(val) for val in center]
            if not skip_validate:
                self.validate()
        except ValueError as err:
            raise LatticeError("Values in center must be of type float.") from err

    def apply(self, domain, geometry, transform=None, **kwargs):
        """
        Fill a domain with particles within the lattice restricted by the geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param geometry: Geometry defining the region within the lattice in
            which particles are restricted to.
        :type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry

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
        super().__init__(center, skip_validate=True)

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

    def apply(self, domain, geometry, transform=None, **kwargs):
        """
        Fill a domain with particles within the cartesian lattice restricted by the geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param geometry: Geometry defining the region within the lattice in
            which particles are restricted to.
        :type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry

        :param transform: Transformation function applied to each particle.
        :type transform: function

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.
        """
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if not (isinstance(geometry, (Geometry, CombinatoryGeometry)) or \
            type(geometry).__name__ in ('Geometry', 'CombinatoryGeometry')):
            raise LatticeError(
                "geometry must be of type spatialpy.Geometry or spatialpy.CombinatoryGeometry."
            )
        if transform is not None and not isinstance(transform, 'function'):
            raise LatticeError("transform must be a function.")

        count = 0
        for x in numpy.arange(self.xmin, self.xmax + self.deltax, self.deltax):
            for y in numpy.arange(self.ymin, self.ymax + self.deltay, self.deltay):
                for z in numpy.arange(self.zmin, self.zmax + self.deltaz, self.deltaz):
                    if geometry.inside((x, y, z), False):
                        if transform is None:
                            point = [x, y, z]
                        else:
                            point = transform([x, y, z])
                        domain.add_point(point, **kwargs)
                        count += 1
        return count

    def validate(self):
        """
        Validate the cartesian lattice dependencies.
        """
        super().validate()

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
        if not isinstance(self.deltax, (int, float)):
            raise LatticeError("deltax must be of type float.")
        if self.deltax <= 0:
            raise LatticeError("deltax must be greater than 0.")
        if not isinstance(self.deltay, (int, float)):
            raise LatticeError("deltay must be of type float.")
        if self.deltay <= 0:
            raise LatticeError("deltay must be greater than 0.")
        if not isinstance(self.deltaz, (int, float)):
            raise LatticeError("deltaz must be of type float.")
        if self.deltaz <= 0:
            raise LatticeError("deltaz must be greater than 0.")

class SphericalLattice(Lattice):
    """
    Spherical lattice class provides a method for creating parts of the spatial
        domain within a spherical coordinate system.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :param radius: Distance between the center and the surface.
    :type radius: float

    :param deltas: Distance between two particle on the surface.
    :type deltas: float

    :param deltar: Radial distance between two particles.
    :type deltar: float

    :raises LatticeError: if center is not and list, doesn't contain 3 values,
        or any value is not a float.
    """
    def __init__(self, radius, deltas, center=None, deltar=None):
        super().__init__(center, skip_validate=True)

        self.radius = radius

        if deltar is None:
            deltar = deltas

        self.deltas = deltas
        self.deltar = deltar

        self.validate()

    def apply(self, domain, geometry, transform=None, **kwargs):
        """
        Fill a domain with particles within the spherical lattice restricted by the geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param geometry: Geometry defining the region within the lattice in
            which particles are restricted to.
        :type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry

        :param transform: Transformation function applied to each particle.
        :type transform: function

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.
        """
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if not (isinstance(geometry, (Geometry, CombinatoryGeometry)) or \
            type(geometry).__name__ in ('Geometry', 'CombinatoryGeometry')):
            raise LatticeError(
                "geometry must be of type spatialpy.Geometry or spatialpy.CombinatoryGeometry."
            )
        if transform is not None and not isinstance(transform, 'function'):
            raise LatticeError("transform must be a function.")

        count = 0
        radius = self.radius
        while radius >= 0:
            # Calculate the approximate number of particle with the radius
            approx_rc = round((4 * radius ** 2) / ((self.deltas / 2) ** 2))

            # Set constants for the radius
            p_area = 4 * numpy.pi * radius ** 2 / approx_rc
            d_a = numpy.sqrt(p_area)
            m_phi = round(numpy.pi * radius / d_a)
            d_phi = numpy.pi / m_phi
            d_theta = p_area / d_phi
            
            for mphi in range(m_phi):
                phi = numpy.pi * (mphi + 0.5) / m_phi
                m_theta = round(2 * numpy.pi * numpy.sin(phi) / d_phi)
            
                for mtheta in range(m_theta):
                    theta = 2 * numpy.pi * mtheta / m_theta
                    x = radius * numpy.cos(theta) * numpy.sin(phi) + self.center[0]
                    y = radius * numpy.sin(theta) * numpy.sin(phi) + self.center[1]
                    z = radius * numpy.cos(phi) + self.center[2]
                    if geometry.inside((x, y, z), False):
                        if transform is None:
                            point = [x, y, z]
                        else:
                            point = transform([x, y, z])
                        domain.add_point(point, **kwargs)
                        count += 1
            radius -= self.deltar
        return count

    def validate(self):
        """
        Validate the spherical lattice dependencies.
        """
        super().validate()

        if not isinstance(self.radius, (int, float)):
            raise LatticeError("radius must be of type float.")
        if self.radius <= 0:
            raise LatticeError("radius must be greater than 0.")
        if not isinstance(self.deltas, (int, float)):
            raise LatticeError("deltas must be of type float.")
        if self.deltas <= 0:
            raise LatticeError("deltas must be greater than 0.")
        if not isinstance(self.deltar, (int, float)):
            raise LatticeError("deltar must be of type float.")
        if self.deltar <= 0:
            raise LatticeError("deltar must be greater than 0.")
