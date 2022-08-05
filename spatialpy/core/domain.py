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
from itertools import combinations
from collections import OrderedDict

import numpy
from plotly.offline import init_notebook_mode, iplot
from scipy.spatial import KDTree

from spatialpy.core.geometry import GeometryAll
from spatialpy.core.lattice import Lattice
from spatialpy.core.transformation import Transformation
from spatialpy.core.visualization import Visualization
from spatialpy.core.spatialpyerror import DomainError

class Domain():
    """
    Domain class for SpatialPy.  A domain defines points and attributes of a regional space for simulation.

    :param numpoints: Total number of spatial domain points
    :type numpoints: int

    :param xlim: Range of domain along x-axis
    :type xlim: float(2)

    :param ylim: Range of domain along y-axis
    :type ylim: float(2)

    :param zlim: Range of domain along z-axis
    :type zlim: float(2)

    :param rho: Background density for the system
    :type rho: float

    :param c0: Speed of sound for the system
    :type c0: float

    :param P0: Background pressure for the system
    :type P0: float

    :param gravity: Acceleration of gravity for the system.
    :type gravity: float[3]
    """
    def __init__(self, numpoints, xlim, ylim, zlim, rho0=1.0, c0=10, P0=None, gravity=None):
        self.vertices = numpy.zeros((numpoints, 3), dtype=float)
        self.triangles = None
        self.tetrahedrons = None

        self.actions = []
        self.on_boundary = None
        self.domain_size = None
        self.tetrahedron_vol = None
        self.dimensions = None

        self.vol = numpy.zeros((numpoints), dtype=float)
        self.mass = numpy.zeros((numpoints), dtype=float)
        self.type_id = numpy.array(["type_UnAssigned"] * numpoints, dtype=object)
        self.nu = numpy.zeros((numpoints), dtype=float)
        self.c = numpy.zeros((numpoints), dtype=float)
        self.rho = numpy.zeros((numpoints), dtype=float)
        self.fixed = numpy.zeros((numpoints), dtype=bool)
        self.listOfTypeIDs = []
        self.typeNdxMapping = OrderedDict({"type_UnAssigned": 0})
        self.typeNameMapping = None

        self.rho0 = rho0
        self.c0 = c0
        if P0 is not None:
            self.P0 = P0
        else:
            self.P0 = rho0 * c0**2 # assume gamma of 1, full equation is rho0*c0**2/gamma
        self.gravity = gravity

        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim

    def __str__(self):
        pad = "    "
        domain_strs = ["Domain Attributes", "", f"{pad}On Boundary: {self.on_boundary}",
                     f"{pad}Domain Size: {self.domain_size}", f"{pad}RHO_0: {self.rho0}", f"{pad}C_0: {self.c0}",
                     f"{pad}P_0: {self.P0}", f"{pad}Gravity: {self.gravity}", f"{pad}X Limit: {self.xlim}",
                     f"{pad}Y Limit: {self.ylim}", f"{pad}Z Limit: {self.zlim}"]

        print("\n".join(domain_strs))

        self._ipython_display_()

        return ""

    def __set_particle_properties(self, ndx, type_id=None, vol=None,
                                  mass=None, nu=None, rho=None, c=None, fixed=False):
        try:
            if type_id is not None
                self.type_id[ndx] = str(type_id)
            if vol is not None:
                self.vol[ndx] = float(vol)
            if mass is not None:
                self.mass[ndx] = float(mass)
            if rho is not None:
                self.rho[ndx] = float(rho)
            if nu is not None:
                self.nu[ndx] = float(nu)
            if c is not None:
                self.c[ndx] = float(c)
            self.fixed[ndx] = bool(fixed)
        except ValueError as err:
            raise DomainError(f"Failed to set all properties. Reaseon given: {err}") from err

    def _ipython_display_(self, use_matplotlib=False):
        self.plot_types(width="auto", height="auto", use_matplotlib=use_matplotlib)

    def _get_type_name_mapping(self):
        self.typeNameMapping = OrderedDict()
        for name, ndx in self.typeNdxMapping.items():
            self.typeNameMapping[ndx] = name

    def apply_actions(self, start=0, end=None, preserve_actions=False):
        """
        Apply domain actions in order from start to end.

        :param start: Starting index for actions (inclusive).
        :type start: int

        :param end: Ending index for actions (exclusive).
        :type end: int

        :param preserve_actions: If False, remove the action after its applied.
        :type preserve_actions: bool

        :raises DomainError: If the action is not supported.
        """
        if end is None:
            end = len(self.actions)
        
        count = end - start
        while count > 0:
            action = self.actions[start]
            if action['type'] == "fill":
                self.apply_fill_action(action)
            elif action['type'] == "set":
                self.apply_set_action(action)
            elif action['type'] == "remove":
                self.apply_remove_action(action)
            else:
                raise DomainError(f"Action of type {action['type']} is not currently supported.")

            if preserve_actions:
                start += 1
            else:
                self.actions.pop(start)
            count -= 1

    def apply_fill_action(self, action):
        """
        Add particles within a region defined by the actions lattice and geometry to the domain.
        Particles will have attributes defined in the actions props.

        :param action: Fill action containing a lattice, geometry, and particle properties.
            Example: {'type': 'fill', 'lattice': lattice_obj, 'geometry':geometry_obj, 'props':{}}
        :type action: dict

        :raises DomainError: If action is not a fill action or is missing a lattice.
        """
        if action['geometry'] is None:
            action['geometry'] = GeometryAll()

        self.validate_action(action, "fill")

        action['lattice'].apply(self, action['geometry'], **action['props'])

    def apply_remove_action(self, action):
        """
        Remove particles within a domain region defined by the actions geometry.

        :param action: Remove action containing a geometry.
            Example: {'type': 'remove', 'geometry':geometry_obj}
        :type action: dict

        :raises DomainError: If action is not a set action, the actions geometry is an invalid type,
            the action is missing props, or type_id in props contains an invalid character or is an int < 0.
        """
        if action['geometry'] is None:
            action['geometry'] = GeometryAll()

        self.validate_action(action, "remove")

        # remove the particles that fall within the defined region
        on_boundary = self.find_boundary_points(update=True)
        for v_ndx in range(self.get_num_voxels()):
            if action['geometry'].inside(self.coordinates()[v_ndx, :], on_boundary[v_ndx]):
                self.type_id = numpy.delete(self.type_id, v_ndx, 0)
                self.vol = numpy.delete(self.vol, v_ndx)
                self.mass = numpy.delete(self.mass, v_ndx)
                self.rho = numpy.delete(self.rho, v_ndx)
                self.nu = numpy.delete(self.nu, v_ndx)
                self.c = numpy.delete(self.c, v_ndx)
                self.fixed = numpy.delete(self.fixed, v_ndx)
                self.on_boundary = numpy.delete(self.on_boundary, v_ndx)

    def apply_set_action(self, action):
        """
        Set properties of particles within a domain region defined by the actions geometry.
        Particles will have attributes defined in the actions props.

        :param action: Set action containing a geometry and particle properties.
            Example: {'type': 'set', 'geometry':geometry_obj , 'props':{}}
        :type action: dict

        :raises DomainError: If action is not a set action, the actions geometry is an invalid type,
            the action is missing props, or type_id in props contains an invalid character or is an int < 0.
        """
        if action['geometry'] is None:
            action['geometry'] = GeometryAll()

        self.validate_action(action, "set")

        if "type_id" in action['props']:
            if isinstance(action['props']['type_id'], int) and action['props']['type_id'] <= 0:
                raise DomainError("Type_id must be a non-zero positive integer or a string.")
            action['props']['type_id'] = f"type_{action['props']['type_id']}"
            for char in action['props']['type_id']:
                if (char in string.punctuation and char != "_") or char == " ":
                    raise DomainError(f"Type_id cannot contain '{char}'")
            if action['props']['type_id'] not in self.typeNdxMapping:
                if "UnAssigned" in action['props']['type_id']:
                    self.typeNdxMapping[action['props']['type_id']] = 0
                else:
                    self.typeNdxMapping[action['props']['type_id']] = len(self.typeNdxMapping)
        
        # apply the properties to all points that fall within the defined region
        on_boundary = self.find_boundary_points(update=True)
        for v_ndx in range(self.get_num_voxels()):
            if action['geometry'].inside(self.coordinates()[v_ndx, :], on_boundary[v_ndx]):
                self.__set_particle_properties(ndx, **action['props'])

    def validate_action(self, action, coverage):
        """
        Validate a domain action.

        :param action: Domain action to be validated.
        :type action: dict

        :param coverage: Scope of the validation.  Accepted values: 'fill', 'set', 'remove'.
        :type coverage: str

        :raises DoaminError: If one of the following conditions are met: The action is an invalid type.
            The action's geometry, lattice or props are not a valid type.
        """
        if not isinstance(action, dict):
            raise DomainError("Actions must be of type dict.")

        if coverage in ("set", "remove"):
            g_type = type(action['geometry']).__name__
            if not (isinstance(action['geometry'], (Geometry, Transformation)) or \
                                                    g_type in ("Geometry", "Transformation")):
                raise DomainError(
                    f"A {action['type']} action's geometry must be of type 'Geometry' or 'Transformation' not {g_type}"
                )

        if coverage == "fill":
            if action['type'] != "fill":
                raise DomainError(f"The action's type must be 'fill' not '{action['type']}'.")

            if action['lattice'] is None:
                raise DomainError("Fill actions must have a lattice.")
            l_type = type(action['lattice']).__name__
            if not (isinstance(action['lattice'], (Lattice, Transformation)) or \
                                                    l_type in ("Lattice", "Transformation")):
                raise DomainError(
                    f"A fill action's lattice must be of type 'Lattice' or 'Transformation' not {l_type}"
                )

        if coverage == "set":
            if action['type'] != "set":
                raise DomainError(f"The action's type must be 'set' not '{action['type']}'.")

            if action['props'] is None or action['props'] == {}:
                raise DomainError("Set actions must have a props.")

        if coverage == "remove":
            if action['type'] != "remove":
                raise DomainError(f"The action's type must be 'remove' not '{action['type']}'.")

        if coverage in ("fill", "set"):
            if action['props'] is not None and not isinstance(action['props'], dict):
                raise DomainError(f"An action's kwargs must be of type dict not {type(action['props'])}")

    def get_type_def(self, type_id):
        """
        Get the C++ type definition for the given type.

        :param type_id: The type_id within the domain.
        :type type_id: str

        :returns: The C++ type definition for the type_id.
        :rtype: str

        :raises DomainError: If the type is not defined within the domain.
        """
        type_id = f"type_{type_id}"
        if type_id not in self.typeNdxMapping:
            errmsg = f"Type_id {type_id} could not be found. "
            errmsg += "Use Domain.set_properties to set the type_id for particles."
            raise DomainError(errmsg)
        return type_id

    def compile_prep(self):
        """
        Generate the domains list of type ids and check for invalid type_ids and rho values
        in preperation of compiling the simulation files.

        :raises DomainError: If a type_id is not set or rho=0 for a particle.
        """
        if self.type_id.tolist().count("type_UnAssigned") > 0:
            raise DomainError(f"Particles must be assigned a type_id.")
        if numpy.count_nonzero(self.rho) < len(self.rho):
            raise DomainError(f"Rho must be a positive value.")

        self.listOfTypeIDs = list(self.typeNdxMapping.values())
        self._get_type_name_mapping()

    def add_point(self, point, vol=1, mass=1, type_id=1, nu=0, fixed=False, rho=None, c=10):
        """
        Add a single point particle to the domain space.

        :param point: Spatial coordinate vertices of point to be added
        :type point: float(3)

        :param vol: Default volume of particle to be added
        :type vol: float

        :param mass: Default mass of particle to be added
        :type mass: float

        :param type_id: Particle type ID of particle to be craddedeated
        :type type_id: str | int

        :param nu: Default viscosity of particle to be added
        :type nu: float

        :param fixed: True if particle is spatially fixed, else False
        :type fixed: bool

        :param c: Default artificial speed of sound of particle to be added
        :type c: float

        :param rho: Default density of particle to be added
        :type rho: float

        :raises DomainError: Type_id is 0 or type_id contains an invalid character.
        """
        if vol < 0:
            raise DomainError("Volume must be a positive value.")

        if isinstance(type_id, int) and type_id <= 0:
            raise DomainError("Type_id must be a non-zero positive integer or a string.")
        type_id = f"type_{type_id}"
        for char in type_id:
            if (char in string.punctuation and char != "_") or char == " ":
                raise DomainError(f"Type_id cannot contain {char}")
        if type_id not in self.typeNdxMapping:
            if "UnAssigned" in type_id:
                self.typeNdxMapping[type_id] = 0
            else:
                self.typeNdxMapping[type_id] = len(self.typeNdxMapping)

        if rho is None:
            rho = mass / vol

        self.type_id = numpy.append(self.type_id, type_id)
        self.vol = numpy.append(self.vol, vol)
        self.mass = numpy.append(self.mass, mass)
        self.nu = numpy.append(self.nu, nu)
        self.c = numpy.append(self.c, c)
        self.rho = numpy.append(self.rho, rho)
        self.fixed = numpy.append(self.fixed, fixed)

        self.vertices = numpy.append(self.vertices, [point], axis=0)

    def set_properties(self, geometry_ivar, type_id, vol=None, mass=None, nu=None, rho=None, c=None, fixed=False):
        """
        Add a type definition to the domain. By default, all regions are set to type 0.

        :param geometry_ivar: an instance of a :py:class:`spatialpy.core.geometry.Geometry` subclass. \
                   The 'inside()' method of this object will be used to assign properties to points.
        :type geometry_ivar: spatialpy.core.geometry.Geometry

        :param type_id: The identifier for this type.
        :type type_id: str | int

        :param vol: The volume of each particle in the type.
        :type vol: float

        :param mass: The mass of each particle in the type.
        :type mass: float

        :param rho: The density of each particle in the type.
        :type rho: float

        :param nu: The viscosity of each particle in the type.
        :type nu: float

        :param c: The artificial speed of sound of each particle in the type.
        :type c: float

        :param fixed: Are the particles in this type immobile.
        :type fixed: bool

        :returns: The number of particles that were tagged with this type_id.
        :rtype: int

        :raises DomainError: Type_id is 0 or type_id contains an invalid character.
        """
        if isinstance(type_id, int) and type_id <= 0:
            raise DomainError("Type_id must be a non-zero positive integer or a string.")
        type_id = f"type_{type_id}"
        for char in type_id:
            if (char in string.punctuation and char != "_") or char == " ":
                raise DomainError(f"Type_id cannot contain '{char}'")
        if type_id not in self.typeNdxMapping:
            if "UnAssigned" in type_id:
                self.typeNdxMapping[type_id] = 0
            else:
                self.typeNdxMapping[type_id] = len(self.typeNdxMapping)
        # apply the type to all points, set type for any points that match
        count = 0
        on_boundary = self.find_boundary_points()
        for v_ndx in range(self.get_num_voxels()):
            if geometry_ivar.inside(self.coordinates()[v_ndx, :], on_boundary[v_ndx]):
                self.type_id[v_ndx] = type_id
                if vol is not None:
                    self.vol[v_ndx] = vol
                if mass is not None:
                    self.mass[v_ndx] = mass
                if rho is not None:
                    self.rho[v_ndx] = rho
                if nu is not None:
                    self.nu[v_ndx] = nu
                if c is not None:
                    self.c[v_ndx] = c
                self.fixed[v_ndx] = fixed
                count +=1
        if count == 0:
            from spatialpy.core import log # pylint: disable=import-outside-toplevel
            log.warning("Type with type_id={} has zero particles in it", type_id)
        return count

    def fill_with_particles(self, geometry_ivar, deltax, deltay=None, deltaz=None, xmin=None,
                            xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, **kwargs):
        """
        Fill a geometric shape with particles.

        :param geometry_ivar: an instance of a :py:class:`spatialpy.core.geometry.Geometry` subclass. \
                   The 'inside()' method of this object will be used to create add the particles.
        :type geometry_ivar: spatialpy.core.geometry.Geometry

        :param deltax: Distance between particles on the x-axis.
        :type deltax: float

        :param deltay: Distance between particles on the y-axis (defaults to deltax).
        :type deltay: float

        :param deltaz: Distance between particles on the z-axis (defaults to deltax).
        :type deltaz: float

        :param xmin: Minimum x value of the bounding box (defaults to Domain.xlim[0]).
        :type xmin: float

        :param xmax: Maximum x value of the bounding box (defaults to Domain.xlim[1]).
        :type xmax: float

        :param ymin: Minimum y value of the bounding box (defaults to Domain.ylim[0]).
        :type ymin: float

        :param ymax: Maximum y value of the bounding box (defaults to Domain.ylim[1]).
        :type ymax: float

        :param zmin: Minimum z value of the bounding box (defaults to Domain.zlim[0]).
        :type zmin: float

        :param zmax: Maximum z value of the bounding box (defaults to Domain.zlim[1]).
        :type zmax: float

        :param \**kwargs: Additional keyword arguments passed to :py:meth:`Domain.add_point`.

        :returns: The number of particles that were created within this geometry.
        :rtype: int
        """
        if deltax <= 0:
            raise DomainError("Deltax must be greater than 0.")
        if deltay is None:
            deltay = deltax
        elif deltay <= 0:
            raise DomainError("Deltay must be greater than 0.")
        if deltaz is None:
            deltaz = deltax
        elif deltaz <= 0:
            raise DomainError("Deltaz must be greater than 0.")

        if xmin is None:
            xmin = self.xlim[0]
        if xmax is None:
            xmax = self.xlim[1]
        if ymin is None:
            ymin = self.ylim[0]
        if ymax is None:
            ymax = self.ylim[1]
        if zmin is None:
            zmin = self.zlim[0]
        if zmax is None:
            zmax = self.zlim[1]

        count = 0
        for x in numpy.arange(xmin, xmax + deltax, deltax):
            for y in numpy.arange(ymin, ymax + deltay, deltay):
                for z in numpy.arange(zmin, zmax + deltaz, deltaz):
                    if geometry_ivar.inside((x, y, z), False):
                        self.add_point([x, y, z], **kwargs)
                        count += 1
        return count

    def find_boundary_points(self, update=False):
        """
        Find all vertices that exist on boundary.

        :returns: A numpy array indexed by vertices, True for boundary points, else false.
        :rtype: np.ndarray(dtype=bool)
        """
        if update or self.on_boundary is None:
            self.on_boundary = numpy.zeros((self.get_num_voxels()), dtype=bool)
            # exterior triangles are part of one-and-only-one tetrahedron
            if self.triangles is None or len(self.triangles) == 0 or len(self.tetrahedrons) == 0:
                return self.on_boundary
            triangle_in_tetrahedrons_count = {}
            for i in range(self.tetrahedrons.shape[0]):
                tets = self.tetrahedrons[i, :]
                tets.sort()
                for p in combinations(tets, 3):
                    key = ".".join([str(s) for s in p])
                    if key in triangle_in_tetrahedrons_count:
                        triangle_in_tetrahedrons_count[key] += 1
                    else:
                        triangle_in_tetrahedrons_count[key] = 1
            boundary_points = set({})
            for key, count in triangle_in_tetrahedrons_count.items():
                if count == 1:
                    (a, b, c) = key.split('.')
                    boundary_points.add(int(a))
                    boundary_points.add(int(b))
                    boundary_points.add(int(c))
            for vertex in boundary_points:
                self.on_boundary[vertex] = True
        return self.on_boundary

    def get_domain_size(self):
        """
        Estimate of domain size at each vertex as the average of the
        diameters of the circumradius of the tetrahedrons that vertex
        is a part of.

        :returns: a numpy array containing the mean for each vertex based on all incident cells
        :rtype: numpy.array
        """
        if self.domain_size is None:
            _ = self.get_vol()

            # Compute the circumradius of the cells
            cr = numpy.zeros((self.tetrahedrons.shape[0]), dtype=float)
            for i, _ in enumerate(cr):
                t_vtx = self.tetrahedrons[i, :]
                # https://en.wikipedia.org/wiki/Tetrahedron#Circumradius
                a = self.distance_between_2_vertices(t_vtx[0], t_vtx[1])
                A = self.distance_between_2_vertices(t_vtx[2], t_vtx[3])
                b = self.distance_between_2_vertices(t_vtx[0], t_vtx[2])
                B = self.distance_between_2_vertices(t_vtx[1], t_vtx[3])
                c = self.distance_between_2_vertices(t_vtx[0], t_vtx[3])
                C = self.distance_between_2_vertices(t_vtx[1], t_vtx[2])
                R = numpy.sqrt((a*A+b*B+c*C)*(a*A+b*B-c*C)*(a*A-b*B+c*C)*(-a*A+b*B+c*C)) / (24*self.tetrahedron_vol[i])
                cr[i] = R

            # Compute the mean for each vertex based on all incident cells
            self.domain_size = numpy.zeros((self.vertices.shape[0]), dtype=float)
            count = numpy.zeros((self.vertices.shape[0]), dtype=float)
            for tndx in range(self.tetrahedrons.shape[0]):
                for vndx in self.tetrahedrons[tndx, :]:
                    self.domain_size[vndx] += cr[tndx]
                    count[vndx] += 1
            for vndx, _ in enumerate(self.domain_size):
                self.domain_size[vndx] = self.domain_size[vndx] / count[vndx]

        return self.domain_size

    def distance_between_2_vertices(self, start, end):
        """
        Get distance between 2 domain vertices.

        :param start: Starting point
        :type start: float(3)

        :param end: Ending point
        :type end: float(2)

        :returns: a distance measurement between start and end point
        :rtype: float
        """
        return numpy.linalg.norm(self.vertices[start, :] - self.vertices[end, :])

    def closest_vertex(self, point):
        """
        Find the nearest vertex of a given point in the domain.

        :param point: Target source point
        :type point: float(3)

        :returns: The coordinates of the nearest vertex to the source point.
        :rtype: float(3)
        """
        min_dist = None
        min_vtx = None
        for i in range(self.vertices.shape[0]):
            d = numpy.linalg.norm(self.vertices[i, :] - point)
            if d > 0 and (min_dist is None or d < min_dist):
                min_dist = d
                min_vtx = i
        return min_vtx

    def coordinates(self):
        """
        Get coordinates within domain.

        :returns: Spatial coordinate vertices of points.
        :rtype: numpy.array
        """
        return self.vertices

    def get_num_voxels(self):
        """
        Get number of voxels in domain.

        :returns: Number of voxels in the domain.
        :rtype: int
        """
        return self.vertices.shape[0]

    def find_h(self):
        """
        Find h value of system.  This value is based off of \
        the particle which has the greatest distance to \
        its nearest neighbor.

        :return: Greatest distance between particle and nearest neighbor.
        :rtype: float
        """
        kdtree = KDTree(self.vertices)
        # Detect nearest neighbor distances for all points
        # since each searched point is already included in
        # the tree, we search 2 nearest neighbors, since
        # the first is just the point itself
        distances, _ = kdtree.query(self.vertices, 2)
        # We only need the distances to the second (non-self) neighbor.
        max_dist = max(distances[:, 1])
        h = 2.2 * max_dist
        return h

    def get_bounding_box(self):
        """
        Get the bounding box of the entire domain.

        :returns: Limits of the bounding box.
        :rtype: tuple(float(2), float(2), float(2))
        """
        xlim = (min(self.vertices[:, 0]), max(self.vertices[:, 0]))
        ylim = (min(self.vertices[:, 1]), max(self.vertices[:, 1]))
        zlim = (min(self.vertices[:, 2]), max(self.vertices[:, 2]))
        return xlim, ylim, zlim

    def get_vol(self):
        """
        Get the total volume of the domain.

        :returns: Total volume of the system.
        :rtype: float
        """
        if self.vol is None:
            self.calculate_vol()
        return self.vol

    def calculate_vol(self):
        """
        Calculate the total volume of the domain.
        """
        self.vol = numpy.zeros((self.vertices.shape[0]), dtype=float)
        self.tetrahedron_vol = numpy.zeros((self.tetrahedrons.shape[0]), dtype=float)
        for t_ndx in range(self.tetrahedrons.shape[0]):
            v1, v2, v3, v4 = self.tetrahedrons[t_ndx]
            a = self.vertices[v1, :]
            b = self.vertices[v2, :]
            c = self.vertices[v3, :]
            d = self.vertices[v4, :]
            #https://en.wikipedia.org/wiki/Tetrahedron#Volume
            t_vol = numpy.abs(numpy.dot((a - d), numpy.cross((b - d), (c - d))) / 6)
            self.tetrahedron_vol[t_ndx] = t_vol
            self.vol[v1] += t_vol / 4
            self.vol[v2] += t_vol / 4
            self.vol[v3] += t_vol / 4
            self.vol[v4] += t_vol / 4

    def plot_types(self, width=None, height=None, colormap=None, size=None, title=None,
                   included_types_list=None, use_matplotlib=False, return_plotly_figure=False):
        '''
        Plots the domain using plotly. Can only be viewed in a Jupyter Notebook.

        :param width: Width in pixels of output plot box or for matplotlib inches of output plot box. \
        Default=500
        :type width: int

        :param height: Height in pixels of output plot box or for matplotlib inches of output plot box. \
        Default=500
        :type height: int

        :param colormap: colormap to use.  Plotly specification. **valid values:** \
        "Plotly3","Jet","Blues","YlOrRd", "PuRd","BuGn","YlOrBr","PuBuGn","BuPu",\
        "YlGnBu", "PuBu","GnBu","YlGn","Greens","Reds", "Greys","RdPu","OrRd","Purples","Oranges".
        :type colormap: str

        :param size: Size in pixels of the particle
        :type size: int

        :param title: The title of the graph
        :type title: str

        :param included_types_list: A list of ints describing which types to include. By default displays all types.
        :type included_types_list: list

        :param return_plotly_figure: Whether or not to return a figure dictionary of data(graph object traces) \
        and layout options which may be edited by the user.
        :type return_plotly_figure: bool

        :param use_matplotlib: Whether or not to plot the proprties results using matplotlib.
        :type use_matplotlib: bool

        :returns: Plotly figure of domain types if, use_matplotlib=False and return_plotly_figure=True
        :rtype: None or dict
        '''
        if len(self.vertices) == 0:
            raise DomainError("The domain does not contain particles.")

        from spatialpy.core.result import _plotly_iterate # pylint: disable=import-outside-toplevel

        if not use_matplotlib:
            if width in (None, "auto"):
                width = None if width == "auto" else 500
            if height in (None, "auto"):
                height = None if height == "auto" else 500

        if not numpy.count_nonzero(self.vertices[:, 1]):
            self.dimensions = 1
        elif not numpy.count_nonzero(self.vertices[:, 2]):
            self.dimensions = 2
        else:
            self.dimensions = 3

        self._get_type_name_mapping()

        types = {}
        # Normalize volumes to [0, 1]
        ptp = numpy.ptp(self.vol)
        if ptp == 0:
            vols = numpy.array([0.5] * len(self.vol))
        else:
            vols = (self.vol - numpy.min(self.vol)) / ptp
        for i, type_id in enumerate(self.type_id):
            name = type_id[5:]
            if included_types_list is None or name in included_types_list:
                if name in types:
                    types[name]['points'].append(self.vertices[i])
                    types[name]['data'].append(self.typeNdxMapping[type_id])
                    types[name]['size_scale'] = numpy.append(types[name]['size_scale'], vols[i])
                else:
                    types[name] = {
                        "points": [self.vertices[i]],
                        "data": [self.typeNdxMapping[type_id]],
                        "size_scale": numpy.array([vols[i]])
                    }

        if use_matplotlib:
            if not isinstance(use_matplotlib, dict):
                use_matplotlib = {}
            use_matplotlib['limits'] = (
                (self.xlim[0] - 0.25, self.xlim[1] + 0.25), (self.ylim[0] - 0.25, self.ylim[1] + 0.25)
            )

            # Support for width, height, and title args
            if width not in (None, "auto") and height not in (None, "auto"):
                # TODO: Deprecation warning for width and height
                plot_args = {"figsize": (width, height)}

                if "plot_args" in use_matplotlib:
                    for name, val in use_matplotlib['plot_args'].items():
                        plot_args[name] = val
                use_matplotlib['plot_args'] = plot_args

            base_group_args = {}
            if colormap is not None:
                base_group_args['cmap'] = colormap
                base_group_args['vmin'] = 1 # minimum number of defined types
                base_group_args['vmax'] = len(self.typeNdxMapping) # number of defined types
            if size is not None:
                base_group_args['s'] = size

            if "scatter_args" not in use_matplotlib:
                use_matplotlib['scatter_args'] = {}
            for type_id in self.typeNdxMapping.keys():
                type_id = type_id[5:]
                group_args = base_group_args.copy()
                if type_id in use_matplotlib['scatter_args']:
                    for name, val in use_matplotlib['scatter_args'][type_id].items():
                        group_args[name] = val
                use_matplotlib['scatter_args'][type_id] = group_args

            if title is not None:
                use_matplotlib['title'] = title

            vis_obj = Visualization(data=types)
            vis_obj.plot_scatter(**use_matplotlib)
            return

        if size is None:
            size = 5

        is_2d = self.dimensions == 2

        trace_list = _plotly_iterate(types, size=size, property_name="type",
                                     colormap=colormap, is_2d=is_2d)

        scene = {
            "aspectmode": 'data',
        }
        layout = {"width": width, "height": width, "scene":scene,
                  "xaxis":{"range":self.xlim}, "yaxis":{"range":self.ylim}
                 }

        if title is not None:
            layout["title"] = title

        fig = {"data":trace_list, "layout":layout}

        if return_plotly_figure:
            return fig
        init_notebook_mode(connected=True)
        iplot(fig)
        return

    @classmethod
    def read_xml_mesh(cls, filename):
        """
        Read a FEniCS/dolfin style XML mesh file

        :param filename: Name of file to read.
        :type filename: str

        :returns: SpatialPy Domain object created from xml mesh.
        :rtype: spatialpy.core.domain.Domain
        """
        root = ET.parse(filename).getroot()
        if not root.tag == 'dolfin':
            raise DomainError("Not a FEniCS/dolfin xml mesh.")
        mesh = root[0]
        if mesh.tag != 'mesh' or \
           mesh.attrib['celltype'] != 'tetrahedron' or \
           mesh.attrib['dim'] != '3':
            raise DomainError("XML mesh format error")
        vertices = mesh[0]
        cells = mesh[1]
        #vertices
        mesh_vertices = numpy.zeros(( len(vertices), 3), dtype=float)
        for vertex in vertices:
            mesh_vertices[int(vertex.attrib['index']), 0] = float(vertex.attrib['x'])
            mesh_vertices[int(vertex.attrib['index']), 1] = float(vertex.attrib['y'])
            mesh_vertices[int(vertex.attrib['index']), 2] = float(vertex.attrib['z'])

        # create domain object
        xlim = (min(mesh_vertices[:, 0]), max(mesh_vertices[:, 0]) )
        ylim = (min(mesh_vertices[:, 1]), max(mesh_vertices[:, 1]) )
        zlim = (min(mesh_vertices[:, 2]), max(mesh_vertices[:, 2]) )
        obj = Domain(len(vertices), xlim, ylim, zlim)
        obj.vertices = mesh_vertices

        #tetrahedrons
        obj.tetrahedrons = numpy.zeros((len(cells), 4), dtype=int)
        for cell in cells:
            obj.tetrahedrons[int(cell.attrib['index']), 0] = int(cell.attrib['v0'])
            obj.tetrahedrons[int(cell.attrib['index']), 1] = int(cell.attrib['v1'])
            obj.tetrahedrons[int(cell.attrib['index']), 2] = int(cell.attrib['v2'])
            obj.tetrahedrons[int(cell.attrib['index']), 3] = int(cell.attrib['v3'])
        # volume
        obj.calculate_vol()
        if not numpy.count_nonzero(obj.vol):
            raise DomainError("Paritcles cannot have 0 volume")
        # set Mass equal to the volume
        obj.mass = obj.vol
        # Calculate density
        obj.rho = obj.mass / obj.vol
        # return model ref
        return obj

    @classmethod
    def import_meshio_object(cls, mesh_obj):
        """
        Import a python meshio mesh object.

        :param mesh_obj: MeshIO object to import
        :type mesh_obj: meshio.Mesh

        :returns: SpatialPy Domain object created from the meshio object
        :rtype: spatialpy.core.domain.Domain
        """
        # create domain object
        xlim = (min(mesh_obj.points[:, 0]), max(mesh_obj.points[:, 0]))
        ylim = (min(mesh_obj.points[:, 1]), max(mesh_obj.points[:, 1]))
        zlim = (min(mesh_obj.points[:, 2]), max(mesh_obj.points[:, 2]))
        obj = Domain(len(mesh_obj.points), xlim, ylim, zlim)
        #vertices
        obj.vertices = mesh_obj.points
        # triangles
        triangles = list(filter(lambda cell: cell.type == "triangle", mesh_obj.cells))
        if triangles:
            obj.triangles = triangles[0].data
        #tetrahedrons
        tetras = list(filter(lambda cell: cell.type == "tetra", mesh_obj.cells))
        if tetras:
            obj.tetrahedrons = tetras[0].data
        # volume
        obj.calculate_vol()
        if not numpy.count_nonzero(obj.vol):
            raise DomainError("Paritcles cannot have 0 volume")
        # set Mass equal to the volume
        obj.mass = obj.vol
        # Calculate density
        obj.rho = obj.mass / obj.vol
        # return model ref
        return obj

    @classmethod
    def read_msh_file(cls, filename):
        """
        Read a Gmsh style .msh file

        :param filename: Filename of gmsh file
        :type filename: str

        :returns: SpatialPy Domain object created from the mesh file.
        :rtype: spatialpy.core.domain.Domain
        """
        try:
            import meshio # pylint: disable=import-outside-toplevel
        except ImportError as err:
            raise DomainError("The python package 'meshio' is not installed.") from err

        return cls.import_meshio_object(meshio.read(filename))

    def read_stochss_subdomain_file(self, filename, type_ids=None):
        """
        Read a .txt file that conains the StochSS v1.x spatial subdomain descriptions.

        :param filename: StochSS v1.x subdomain description filename.
        :type filename: str

        :param type_ids: Mapping of type indecies to type names.
        :type type_ids: dict{str:str}

        :raises DomainError: Domain file could not be read or type_id contains an invalid character.
        """
        with open(filename,'r', encoding="utf-8") as file_obj:
            for lnum, line in enumerate(file_obj):
                try:
                    (ndx, type_id) = line.rstrip().split(',')

                    type_id = f"type_{type_id if type_ids is None else type_ids[type_id]}"
                    for char in type_id:
                        if (char in string.punctuation and char != "_") or char == " ":
                            raise DomainError(f"Type_id cannot contain {char}")
                    if type_id not in self.typeNdxMapping:
                        if "UnAssigned" in type_id:
                            self.typeNdxMapping[type_id] = 0
                        else:
                            self.typeNdxMapping[type_id] = len(self.typeNdxMapping)

                    self.type_id[int(ndx)] = type_id

                except ValueError as err:
                    raise DomainError(f"Could not read in subdomain file, error on line {lnum}: {line}") from err

    @classmethod
    def read_stochss_domain(cls, filename):
        """
        Read a StochSS Domain (.domn) file or pull a StochSS Domain from a StochSS Spatial Model (.smdl) file.

        :param filename: Name of file to read.
        :type filename: str

        :returns: SpatialPy Domain object created from StochSS domain.
        :rtype: spatialpy.core.domain.Domain
        """
        try:
            with open(filename, "r", encoding="utf-8") as domain_file:
                domain = json.load(domain_file)
                if "domain" in domain.keys():
                    domain = domain['domain']

            obj = Domain(0, tuple(domain['x_lim']), tuple(domain['y_lim']), tuple(domain['z_lim']),
                        rho0=domain['rho_0'], c0=domain['c_0'], P0=domain['p_0'], gravity=domain['gravity'])

            for i, particle in enumerate(domain['particles']):
                try:
                    type_id = list(filter(
                        lambda d_type, t_ndx=particle['type']: d_type['typeID'] == t_ndx, domain['types']
                    ))[0]['name']
                except IndexError:
                    type_id = particle['type']
                if type_id == "Un-Assigned" or type_id == 0:
                    type_id = "UnAssigned"
                # StochSS backward compatability check for rho
                rho = None if "rho" not in particle.keys() else particle['rho']
                # StochSS backward compatability check for c
                c = 0 if "c" not in particle.keys() else particle['c']
                obj.add_point(particle['point'], vol=particle['volume'], mass=particle['mass'],
                              type_id=type_id, nu=particle['nu'], fixed=particle['fixed'], rho=rho, c=c)

            return obj
        except KeyError as err:
            raise DomainError("The file is not a StochSS Domain (.domn) or a StochSS Spatial Model (.smdl).") from err

    @classmethod
    def create_3D_domain(cls, xlim, ylim, zlim, nx, ny, nz, type_id=1, mass=1.0,
                         nu=1.0, rho=None, c=0, fixed=False, **kwargs):
        """
        Create a filled 3D domain

        :param xlim: highest and lowest coordinate in the x dimension
        :type xlim: float(2)

        :param ylim: highest and lowest coordinate in the y dimension
        :type ylim: float(2)

        :param zlim: highest and lowest coordinate in the z dimension
        :type zlim: float(2)

        :param nx: number of particle spacing in the x dimension
        :type nx: int

        :param ny: number of particle spacing in the y dimension
        :type ny: int

        :param nz: number of particle spacing in the z dimension
        :type nz: int

        :param type_id: default type ID of particles to be created. Defaults to 1
        :type type_id: int

        :param mass: default mass of particles to be created. Defaults to 1.0
        :type mass: float

        :param nu: default viscosity of particles to be created. Defaults to 1.0
        :type nu: float

        :param c: default artificial speed of sound of particles to be created. Defaults to 0.0.
        :type c: float

        :param rho: default density of particles to be created.
        :type rho: float

        :param fixed: spatially fixed flag of particles to be created. Defaults to false.
        :type fixed: bool

        :param \**kwargs: Additional keyword arguments passed to :py:class:`Domain`.

        :returns: Uniform 3D SpatialPy Domain object.
        :rtype: spatialpy.core.domain.Domain
        """
        # Create domain object
        numberparticles = nx * ny * nz
        obj = Domain(0, xlim, ylim, zlim, **kwargs)
        # Vertices
        x_list = numpy.linspace(xlim[0], xlim[1], nx)
        y_list = numpy.linspace(ylim[0], ylim[1], ny)
        z_list = numpy.linspace(zlim[0], zlim[1], nz)
        totalvolume = abs(xlim[1] - xlim[0]) * abs(ylim[1] - ylim[0]) * abs(zlim[1] - zlim[0])
        vol = totalvolume / numberparticles
        if vol < 0:
            raise DomainError("Paritcles cannot have 0 volume")
        for x in x_list:
            for y in y_list:
                for z in z_list:
                    obj.add_point([x, y, z], vol=vol, mass=mass, rho=rho,
                                  type_id=type_id, nu=nu, c=c, fixed=fixed)
        return obj

    @classmethod
    def create_2D_domain(cls, xlim, ylim, nx, ny, type_id=1, mass=1.0,
                         nu=1.0, rho=None, c=0, fixed=False, **kwargs):
        """
        Create a filled 2D domain

        :param xlim: highest and lowest coordinate in the x dimension
        :type xlim: float(2)

        :param ylim: highest and lowest coordinate in the y dimension
        :type ylim: float(2)

        :param nx: number of particle spacing in the x dimension
        :type nx: int

        :param ny: number of particle spacing in the y dimension
        :type ny: int

        :param type_id: default type ID of particles to be created. Defaults to 1
        :type type_id: int

        :param mass: default mass of particles to be created. Defaults to 1.0
        :type mass: float

        :param nu: default viscosity of particles to be created. Defaults to 1.0
        :type nu: float

        :param c: default artificial speed of sound of particles to be created. Defaults to 0.0.
        :type c: float

        :param rho: default density of particles to be created.
        :type rho: float

        :param fixed: spatially fixed flag of particles to be created. Defaults to false.
        :type fixed: bool

        :param \**kwargs: Additional keyword arguments passed to :py:class:`Domain`.

        :returns: Uniform 2D SpatialPy Domain object.
        :rtype: spatialpy.core.domain.Domain
        """
        # Create domain object
        numberparticles = nx * ny
        obj = Domain(0, xlim, ylim, (0, 0), **kwargs)
        # Vertices
        x_list = numpy.linspace(xlim[0], xlim[1], nx)
        y_list = numpy.linspace(ylim[0], ylim[1], ny)
        totalvolume = abs(xlim[1] - xlim[0]) * abs(ylim[1] - ylim[0])
        vol = totalvolume / numberparticles
        if vol < 0:
            raise DomainError("Paritcles cannot have 0 volume")
        for x in x_list:
            for y in y_list:
                obj.add_point([x, y, 0], vol=vol, mass=mass, rho=rho,
                              type_id=type_id, nu=nu, c=c, fixed=fixed)
        return obj
