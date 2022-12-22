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
import json
import string

import xml.etree.ElementTree as ET

import numpy

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
        elif isinstance(center, (list, tuple)):
            center = numpy.array(center)
        try:
            self.center = [float(val) for val in center]
            if not skip_validate:
                self.validate()
        except ValueError as err:
            raise LatticeError("Values in center must be of type float.") from err

    def _update_limits(self, domain):
        xlim, ylim, zlim = domain.get_bounding_box()
        if domain.xlim != xlim:
            domain.xlim = xlim
        if domain.ylim != ylim:
            domain.ylim = ylim
        if domain.zlim != zlim:
            domain.zlim = zlim

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

        self.xmin = xmin + self.center[0]
        self.xmax = xmax + self.center[0]
        self.ymin = ymin + self.center[1]
        self.ymax = ymax + self.center[1]
        self.zmin = zmin + self.center[2]
        self.zmax = zmax + self.center[2]

        if deltay is None:
            deltay = deltax
        if deltaz is None:
            deltaz = deltax

        self.deltax = deltax
        self.deltay = deltay
        self.deltaz = deltaz

        self.validate()

    def __add_point(self, domain, geometry, transform, point, count, kwargs):
        if geometry.inside(point, False):
            if transform is not None:
                point = transform(point)
            domain.add_point(point, **kwargs)
            count += 1
        return count

    def __generate_z(self, domain, geometry, transform, x, y, count, kwargs):
        for z in numpy.arange(self.zmin, self.zmax + self.deltaz, self.deltaz):
            count = self.__add_point(domain, geometry, transform, [x, y, z], count, kwargs)
        return count

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
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if not (isinstance(geometry, (Geometry, CombinatoryGeometry)) or \
            type(geometry).__name__ in ('Geometry', 'CombinatoryGeometry')):
            raise LatticeError(
                "geometry must be of type spatialpy.Geometry or spatialpy.CombinatoryGeometry."
            )
        if transform is not None and not callable(transform):
            raise LatticeError("transform must be a function.")

        count = 0
        for x in numpy.arange(self.xmin, self.xmax + self.deltax, self.deltax):
            for y in numpy.arange(self.ymin, self.ymax + self.deltay, self.deltay):
                if self.deltaz == 0:
                    z = self.center[2]
                    count = self.__add_point(domain, geometry, transform, [x, y, z], count, kwargs)
                else:
                    count = self.__generate_z(domain, geometry, transform, x, y, count, kwargs)
        self._update_limits(domain)
        if 'vol' not in kwargs:
            offset = len(domain.vertices) - count
            if self.deltaz > 0:
                vol = 4 / 3 * numpy.pi * (self.deltax / 2) * (self.deltay / 2) * (self.deltaz / 2)
            else:
                vol = numpy.pi * (self.deltax / 2) * (self.deltay / 2)
            for i in range(offset, offset + count):
                domain.vol[i] = vol
                if 'rho' not in kwargs:
                    domain.rho[i] = domain.mass[i] / vol
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
        if self.deltaz < 0:
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
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if not (isinstance(geometry, (Geometry, CombinatoryGeometry)) or \
            type(geometry).__name__ in ('Geometry', 'CombinatoryGeometry')):
            raise LatticeError(
                "geometry must be of type spatialpy.Geometry or spatialpy.CombinatoryGeometry."
            )
        if transform is not None and not callable(transform):
            raise LatticeError("transform must be a function.")

        count = 0
        radius = self.radius
        while radius > 0:
            # Calculate the approximate number of particle with the radius
            approx_rc = int(round((4 * radius ** 2) / ((self.deltas / 2) ** 2)))

            # Set constants for the radius
            p_area = 4 * numpy.pi * radius ** 2 / approx_rc
            d_a = numpy.sqrt(p_area)
            m_phi = int(round(numpy.pi * radius / d_a))
            d_phi = numpy.pi / m_phi
            d_theta = p_area / d_phi

            for mphi in range(m_phi):
                phi = numpy.pi * (mphi + 0.5) / m_phi
                m_theta = int(round(2 * numpy.pi * numpy.sin(phi) / d_phi))

                for mtheta in range(m_theta):
                    theta = 2 * numpy.pi * mtheta / m_theta
                    x = radius * numpy.cos(theta) * numpy.sin(phi)
                    y = radius * numpy.sin(theta) * numpy.sin(phi)
                    z = radius * numpy.cos(phi)
                    if geometry.inside((x, y, z), False):
                        if transform is None:
                            point = [x, y, z]
                        else:
                            point = transform([x, y, z])
                        if not isinstance(point, numpy.ndarray):
                            point = numpy.array(point)
                        domain.add_point(point + self.center, **kwargs)
                        count += 1
            radius -= self.deltar
        if radius == 0 and geometry.inside((0, 0, 0), False):
            point = [0, 0, 0] if transform is None else transform([0, 0, 0])
            if not isinstance(point, numpy.ndarray):
                point = numpy.array(point)
            domain.add_point(point + self.center, **kwargs)
            count += 1
        self._update_limits(domain)
        if 'vol' not in kwargs:
            offset = len(domain.vertices) - count
            vol = 4 / 3 * numpy.pi * (self.deltas / 2)**2 * (self.deltar / 2)
            for i in range(offset, offset + count):
                domain.vol[i] = vol
                if 'rho' not in kwargs:
                    domain.rho[i] = domain.mass[i] / vol
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

class CylindricalLattice(Lattice):
    """
    Cylindrical lattice class provides a method for creating parts of the spatial
        domain within a cylindrical coordinate system.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :param length: Length of the cylindrical lattice.
    :type length: float

    :param radius: Distance between the center and the surface.
    :type radius: float

    :param deltas: Distance between two particle on the surface.
    :type deltas: float

    :param deltar: Radial distance between two particles.
    :type deltar: float

    :raises LatticeError: if center is not and list, doesn't contain 3 values,
        or any value is not a float.
    """
    def __init__(self, radius, length, deltas, center=None, deltar=None):
        super().__init__(center, skip_validate=True)

        self.radius = radius
        self.length = length

        if deltar is None:
            deltar = deltas

        self.deltas = deltas
        self.deltar = deltar

        self.validate()

    def apply(self, domain, geometry, transform=None, **kwargs):
        """
        Fill a domain with particles within the cylindrical lattice restricted by the geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param geometry: Geometry defining the region within the lattice in
            which particles are restricted to.
        :type geometry: spatialpy.Geometry | spatialpy.CombinatoryGeometry

        :param transform: Transformation function applied to each particle.
        :type transform: function

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.
        """
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if not (isinstance(geometry, (Geometry, CombinatoryGeometry)) or \
            type(geometry).__name__ in ('Geometry', 'CombinatoryGeometry')):
            raise LatticeError(
                "geometry must be of type spatialpy.Geometry or spatialpy.CombinatoryGeometry."
            )
        if transform is not None and not callable(transform):
            raise LatticeError("transform must be a function.")

        count = 0
        h_len = self.length / 2
        xmin = -h_len
        xmax = h_len
        radius = self.radius
        while radius > 0:
            # Calculate the approximate number of particle with the radius
            approx_rc = int(round((2 * radius * self.length) / ((self.deltas / 2) ** 2)))

            p_area = 2 * numpy.pi * radius * self.length / approx_rc
            d_a = numpy.sqrt(p_area)
            m_theta = int(round(2 * numpy.pi * radius / d_a))
            d_theta = 2 * numpy.pi / m_theta

            x = xmin
            while x <= xmax:
                for mtheta in range(m_theta):
                    theta = 2 * numpy.pi * (mtheta + 0.5) / m_theta
                    y = radius * numpy.cos(theta)
                    z = radius * numpy.sin(theta)
                    if geometry.inside((x, y, z), False):
                        if transform is None:
                            point = [x, y, z]
                        else:
                            point = transform([x, y, z])
                        if not isinstance(point, numpy.ndarray):
                            point = numpy.array(point)
                        domain.add_point(point + self.center, **kwargs)
                        count += 1
                x += self.deltas
            radius -= self.deltar
        if radius == 0:
            x = xmin
            while x <= xmax:
                if geometry.inside((x, 0, 0), False):
                    point = [x, 0, 0] if transform is None else transform([x, 0, 0])
                    if not isinstance(point, numpy.ndarray):
                        point = numpy.array(point)
                    domain.add_point(point + self.center, **kwargs)
                    count += 1
                x += self.deltas
        self._update_limits(domain)
        if 'vol' not in kwargs:
            offset = len(domain.vertices) - count
            vol = 4 / 3 * numpy.pi * (self.deltas / 2)**2 * (self.deltar / 2)
            for i in range(offset, offset + count):
                domain.vol[i] = vol
                if 'rho' not in kwargs:
                    domain.rho[i] = domain.mass[i] / vol
        return count

    def validate(self):
        """
        Validate the cylindrical lattice dependencies.
        """
        super().validate()

        if not isinstance(self.radius, (int, float)):
            raise LatticeError("radius must be of type float.")
        if self.radius <= 0:
            raise LatticeError("radius must be greater than 0.")
        if not isinstance(self.length, (int, float)):
            raise LatticeError("length must be of type float.")
        if self.length <= 0:
            raise LatticeError("length must be greater than 0.")
        if not isinstance(self.deltas, (int, float)):
            raise LatticeError("deltas must be of type float.")
        if self.deltas <= 0:
            raise LatticeError("deltas must be greater than 0.")
        if not isinstance(self.deltar, (int, float)):
            raise LatticeError("deltar must be of type float.")
        if self.deltar <= 0:
            raise LatticeError("deltar must be greater than 0.")

class XMLMeshLattice(Lattice):
    """
    XML mesh lattice class provides a method for creating parts of the spatial
        domain with a mesh defined by a FEniCS/dolfin style XML mesh file.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :param filename: Name of file to read.
    :type filename: str

    :param subdomain_file: StochSS v1.x subdomain description filename (optional).
    :type subdomain_file: str

    :param type_ids: Mapping of type indices to type names (optional).
    :type type_ids: dict{str:str}
    """
    def __init__(self, filename, center=None, subdomain_file=None, type_ids=None):
        super().__init__(center, skip_validate=True)

        self.filename = filename
        self.subdomain_file = subdomain_file
        self.type_ids = type_ids

        self.validate()

    def __get_types(self):
        type_ids = {}
        with open(self.subdomain_file,'r', encoding="utf-8") as file_obj:
            for lnum, line in enumerate(file_obj):
                try:
                    (ndx, type_id) = line.rstrip().split(',')

                    if self.type_ids is not None:
                        type_id = self.type_ids[type_id]
                    type_ids[int(ndx)] = f"type_{type_id}"
                except ValueError as err:
                    errmsg = f"Could not read in subdomain file, error on line {lnum}: {line}"
                    raise LatticeError(errmsg) from err
        return type_ids

    def apply(self, domain, *args, transform=None, **kwargs):
        """
        Fill a domain with particles within the xml mesh lattice un-restricted by a geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param transform: Transformation function applied to each particle.
        :type transform: function

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.
        """
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if transform is not None and not callable(transform):
            raise LatticeError("transform must be a function.")

        if self.subdomain_file is not None:
            type_ids = self.__get_types()
        else:
            type_ids = None

        root = ET.parse(self.filename).getroot()
        if not root.tag == 'dolfin':
            raise LatticeError(f"{self.filename} is not a FEniCS/dolfin xml mesh.")

        mesh = root[0]
        if mesh.tag != 'mesh' or mesh.attrib['celltype'] != 'tetrahedron' or mesh.attrib['dim'] != '3':
            raise LatticeError("XML mesh format error.")

        vertices = mesh[0]
        cells = mesh[1]
        #vertices
        mesh_vertices = numpy.zeros(( len(vertices), 3), dtype=float)
        for vertex in vertices:
            mesh_vertices[int(vertex.attrib['index']), 0] = float(vertex.attrib['x'])
            mesh_vertices[int(vertex.attrib['index']), 1] = float(vertex.attrib['y'])
            mesh_vertices[int(vertex.attrib['index']), 2] = float(vertex.attrib['z'])

        count = 0
        for i, vertex in enumerate(mesh_vertices):
            x = vertex[0]
            y = vertex[1]
            z = vertex[2]
            if transform is None:
                point = [x, y, z]
            else:
                point = transform([x, y, z])
            if type_ids is not None and i in type_ids:
                kwargs['type_id'] = type_ids[i]
            if not isinstance(point, numpy.ndarray):
                point = numpy.array(point)
            domain.add_point(point + self.center, **kwargs)
            count += 1

        #tetrahedrons
        if domain.tetrahedrons is None:
            domain.tetrahedrons = numpy.zeros((len(cells), 4), dtype=int)
            for cell in cells:
                domain.tetrahedrons[int(cell.attrib['index']), 0] = int(cell.attrib['v0'])
                domain.tetrahedrons[int(cell.attrib['index']), 1] = int(cell.attrib['v1'])
                domain.tetrahedrons[int(cell.attrib['index']), 2] = int(cell.attrib['v2'])
                domain.tetrahedrons[int(cell.attrib['index']), 3] = int(cell.attrib['v3'])
        else:
            for cell in cells:
                tetrahedron = [
                    int(cell.attrib['v0']), int(cell.attrib['v1']),
                    int(cell.attrib['v2']), int(cell.attrib['v3'])
                ]
                domain.tetrahedrons = numpy.append(domain.tetrahedrons, [tetrahedron], axis=0)

        self._update_limits(domain)
        return count

    def validate(self):
        """
        Validate the XML mesh lattice dependencies.
        """
        super().validate()

        if not isinstance(self.filename, str):
            raise LatticeError("filename must be of type str.")
        if self.subdomain_file is not None and not isinstance(self.subdomain_file, str):
            raise LatticeError("subdomain_file must be of type str.")
        if self.type_ids is not None:
            for ndx, name in self.type_ids.items():
                if not isinstance(ndx, str):
                    raise LatticeError("Keys in type_ids must be of type str.")
                if not isinstance(name, str):
                    raise LatticeError("Values in type_ids must be of type str.")
                for char in name:
                    if (char in string.punctuation and char != "_") or char == " ":
                        raise LatticeError(f"Values in type_ids cannot contain '{char}'")

class MeshIOLattice(Lattice):
    """
    meshio lattice class provides a method for creating parts of the spatial
        domain with a mesh defined by a Gmsh style .msh mesh file.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :param filename: Name of file to read.
    :type filename: str

    :param subdomain_file: StochSS v1.x subdomain description filename (optional).
    :type subdomain_file: str

    :param type_ids: Mapping of type indices to type names (optional).
    :type type_ids: dict{str:str}
    """
    def __init__(self, filename=None, center=None, mesh=None, subdomain_file=None, type_ids=None):
        super().__init__(center, skip_validate=True)

        self.filename = filename
        self.mesh = mesh
        self.subdomain_file = subdomain_file
        self.type_ids = type_ids

        self.validate()

    def __get_types(self):
        type_ids = {}
        with open(self.subdomain_file,'r', encoding="utf-8") as file_obj:
            for lnum, line in enumerate(file_obj):
                try:
                    (ndx, type_id) = line.rstrip().split(',')

                    if self.type_ids is not None:
                        type_id = self.type_ids[type_id]
                    type_ids[int(ndx)] = f"type_{type_id}"
                except ValueError as err:
                    errmsg = f"Could not read in subdomain file, error on line {lnum}: {line}"
                    raise LatticeError(errmsg) from err
        return type_ids

    def apply(self, domain, *args, transform=None, **kwargs):
        """
        Fill a domain with particles within the mesh IO lattice un-restricted by a geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param transform: Transformation function applied to each particle.
        :type transform: function

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.
        """
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if transform is not None and not callable(transform):
            raise LatticeError("transform must be a function.")

        if self.mesh is None:
            try:
                import meshio # pylint: disable=import-outside-toplevel
            except ImportError as err:
                raise LatticeError("The python package 'meshio' is not installed.") from err

            mesh = meshio.read(self.filename)
        else:
            mesh = self.mesh

        if self.subdomain_file is not None:
            type_ids = self.__get_types()
        else:
            type_ids = None

        num_points = len(domain.vertices)
        #vertices
        count = 0
        for i, vertex in enumerate(mesh.points):
            x = vertex[0]
            y = vertex[1]
            z = vertex[2]
            if transform is None:
                point = [x, y, z]
            else:
                point = transform([x, y, z])
            if type_ids is not None and i in type_ids:
                kwargs['type_id'] = type_ids[i]
            if not isinstance(point, numpy.ndarray):
                point = numpy.array(point)
            domain.add_point(point + self.center, **kwargs)
            count += 1

        # triangles
        triangles = list(filter(lambda cell: cell.type == "triangle", mesh.cells))
        if triangles:
            if domain.triangles is None:
                domain.triangles = triangles[0].data
            else:
                for triangle in triangles[0].data:
                    triangle = triangle + num_points
                    domain.triangles = numpy.append(domain.triangles, [triangle], axis=0)

        #tetrahedrons
        tetras = list(filter(lambda cell: cell.type == "tetra", mesh.cells))
        if tetras:
            if domain.tetrahedrons is None:
                domain.tetrahedrons = tetras[0].data
            else:
                for tetra in tetras[0].data:
                    tetra = tetra + num_points
                    domain.tetrahedrons = numpy.append(domain.tetrahedrons, [tetra], axis=0)
        self._update_limits(domain)
        return count

    def validate(self):
        """
        Validate the meshio lattice dependencies.
        """
        super().validate()

        if self.filename is None and self.mesh is None:
            raise LatticeError("MeshIOLattice requires a msh filename or meshio object.")
        if self.filename is not None and not isinstance(self.filename, str):
            raise LatticeError("filename must be of type str.")
        if self.subdomain_file is not None and not isinstance(self.subdomain_file, str):
            raise LatticeError("subdomain_file must be of type str.")
        if self.type_ids is not None:
            for ndx, name in self.type_ids.items():
                if not isinstance(ndx, str):
                    raise LatticeError("Keys in type_ids must be of type str.")
                if not isinstance(name, str):
                    raise LatticeError("Values in type_ids must be of type str.")
                for char in name:
                    if (char in string.punctuation and char != "_") or char == " ":
                        raise LatticeError(f"Values in type_ids cannot contain '{char}'")

class StochSSLattice(Lattice):
    """
    stochss lattice class provides a method for creating parts of the spatial domain
        with a domain defined by a stochss style .domn domain file or .smdl model file.

    :param center: The center point of the lattice.
    :type center: float[3] | float(3)

    :param filename: Name of file to read.
    :type filename: str
    """
    def __init__(self, filename, center=None):
        super().__init__(center, skip_validate=True)

        self.filename = filename

        self.validate()

    def apply(self, domain, *args, transform=None, **kwargs):
        """
        Fill a domain with particles within the stochss lattice un-restricted by a geometry.

        :param domain: Domain particles are to be added to.
        :type domain: spatialpy.Domain

        :param transform: Transformation function applied to each particle.
        :type transform: function
        """
        from spatialpy.core.domain import Domain # pylint: disable=import-outside-toplevel
        if not (isinstance(domain, Domain) or type(domain).__name__ == 'Domain'):
            raise LatticeError("domain must be of type spatialpy.Domain.")
        if transform is not None and not callable(transform):
            raise LatticeError("transform must be a function.")

        try:
            with open(self.filename, "r", encoding="utf-8") as domain_file:
                s_domain = json.load(domain_file)
                if "domain" in s_domain.keys():
                    s_domain = s_domain['domain']

            domain.rho0 = s_domain['rho_0']
            domain.c0 = s_domain['c_0']
            domain.P0 = s_domain['p_0']
            domain.gravity = s_domain['gravity']

            type_ids = {}
            for s_type in s_domain['types']:
                type_ids[s_type['typeID']] = s_type['name'].replace('-', '')

            count = 0
            for particle in s_domain['particles']:
                kwargs = {
                    "type_id": type_ids[particle['type']],
                    "vol": particle['volume'],
                    "mass": particle['mass'],
                    "rho": None if "rho" not in particle.keys() else particle['rho'],
                    "nu": particle['nu'],
                    "c": 0 if "c" not in particle.keys() else particle['c'],
                    "fixed": particle['fixed']
                }

                x = particle['point'][0]
                y = particle['point'][1]
                z = particle['point'][2]
                if transform is None:
                    point = [x, y, z]
                else:
                    point = transform([x, y, z])
                if not isinstance(point, numpy.ndarray):
                    point = numpy.array(point)
                domain.add_point(point + self.center, **kwargs)
                count += 1

            self._update_limits(domain)
            return count
        except KeyError as err:
            raise LatticeError("The file is not a StochSS Domain (.domn) or a StochSS Spatial Model (.smdl).") from err

    def validate(self):
        """
        Validate the stochss lattice dependencies.
        """
        super().validate()

        if not isinstance(self.filename, str):
            raise LatticeError("filename must be of type str.")
