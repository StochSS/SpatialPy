'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

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

import json

import numpy
from scipy.spatial import KDTree

from spatialpy.core.spatialpyError import DomainError

class Domain():
    """ Domain class for SpatialPy.  A domain defines points and attributes of a regional space for simulation.

        :param numpoints: Total number of spatial domain points
        :type numpoints: int
        
        :param xlim: Range of domain along x-axis
        :type xmax: tuple(float)

        :param ylim: Range of domain along y-axis
        :type ymin: tuple(float)
        
        :param zlim: Range of domain along z-axis
        :type ymax: tuple(float)

        :param rho: Background density for the system
        :type zmin: float
        
        :param c0: Speed of sound for the system
        :type zmax: float

        :param P0: Background pressure for the system
        :type typeid: float

        :param gravity: Acceleration of gravity for the system.
        :type species: float
    """


    def __init__(self, numpoints, xlim, ylim, zlim, rho0=1.0, c0=10, P0=10, gravity=None):
        self.vertices = numpy.zeros((numpoints, 3), dtype=float)
        self.triangles = None
        self.tetrahedrons = None

        self.on_boundary = None
        self.domain_size = None
        self.tetrahedron_vol = None
        self.dimensions = None

        self.vol = numpy.zeros((numpoints), dtype=float)
        self.mass = numpy.zeros((numpoints), dtype=float)
        self.type = numpy.zeros((numpoints), dtype=int)
        self.nu = numpy.zeros((numpoints), dtype=float)
        self.c = numpy.zeros((numpoints), dtype=float)
        self.rho = numpy.zeros((numpoints), dtype=float)
        self.fixed = numpy.zeros((numpoints), dtype=bool)

        self.rho0 = rho0
        self.c0 = c0
        self.P0 = P0
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
        domain_strs.extend(["", "Paritcles", ""])
        for i, vertex in enumerate(self.vertices):
            v_str = f"{pad}{i+1}: {vertex}\n{pad}   Volume:{self.vol[i]}, Mass: {self.mass[i]}, "
            v_str += f"Type: {self.type[i]}, Viscosity: {self.nu[i]}, Density: {self.rho[i]}, "
            v_str += f"Artificial Speed of Sound: {self.c[i]}, Fixed: {self.fixed[i]}"
            v_str += f"Type: {self.type[i]}, Viscosity: {self.nu[i]}, Density: {self.rho[i]}, Fixed: {self.fixed[i]}"
            domain_strs.append(v_str)
        if self.triangles is not None:
            domain_strs.extend(["", "Triangles", ""])
            for i, triangle in enumerate(self.triangles):
                domain_strs.append(f"{pad}{i+1}: {triangle}")
        if self.tetrahedrons is not None:
            domain_strs.extend(["", "Tetrahedrons", ""])
            for i, tetrahedron in enumerate(self.tetrahedrons):
                domain_strs.append(f"{pad}{i+1}: {tetrahedron}, Volume: {self.tetrahedron_vol[i]}")

        return "\n".join(domain_strs)

    def _ipython_display_(self, use_matplotlib=False):
        self.plot_types(width="auto", height="auto", use_matplotlib=use_matplotlib)

    def add_point(self, x, vol, mass, type, nu, fixed, rho=None, c=0):
        """ Add a single point particle to the domain space.

            :param x: Spatial coordinate vertices of point to be added
            :type x: tuple(float, float, float) or tuple(float, float)

            :param vol: Default volume of particle to be added
            :type vol: float

            :param mass: Default mass of particle to be added
            :type mass: float

            :param type: Particle typeid of particle to be created
            :type type: int

            :param nu: Default viscosity of particle to be created
            :type nu: float

            :param fixed: True if particle is spatially fixed, else False
            :type fixed: bool

            :param c: Default artificial speed of sound of particle to be created
            :type c: float

            :param rho: Default density of particle to be created
            :type rho: float
        """

        if vol < 0:
            raise DomainError("Volume must be a positive value.")

        if rho is None:
            rho = mass/vol

        self.vol = numpy.append(self.vol, vol)
        self.mass = numpy.append(self.mass, mass)
        self.type = numpy.append(self.type, type)
        self.nu = numpy.append(self.nu, nu)
        self.c = numpy.append(self.c, c)
        self.rho = numpy.append(self.rho, rho)
        self.fixed = numpy.append(self.fixed, fixed)

        self.vertices = numpy.append(self.vertices, [x], axis=0)

    def find_boundary_points(self):
        if self.on_boundary is None:
            self.on_boundary = numpy.zeros((self.get_num_voxels()), dtype=bool)
            # exterior triangles are part of one-and-only-one tetrahedron
            if self.triangles is None or len(self.triangles) == 0 or len(self.tetrahedrons) == 0:
                return self.on_boundary
            from itertools import combinations
            triangle_in_tetrahedrons_count = {}
            for i in range(self.tetrahedrons.shape[0]):
                tets = self.tetrahedrons[i,:]
                tets.sort()
                for p in combinations(tets,3):
                    key = ".".join([str(s) for s in p ])
                    #print(key)
                    if key in triangle_in_tetrahedrons_count:
                        triangle_in_tetrahedrons_count[key]+=1
                    else:
                        triangle_in_tetrahedrons_count[key]=1
            boundary_points = set({})
            for key in triangle_in_tetrahedrons_count:
                #print(key+" "+str(triangle_in_tetrahedrons_count[key]))
                if triangle_in_tetrahedrons_count[key]==1:
                    (a,b,c) = key.split('.')
                    boundary_points.add(int(a))
                    boundary_points.add(int(b))
                    boundary_points.add(int(c))
            for v in boundary_points:
                self.on_boundary[v] = True
        return self.on_boundary


    def get_domain_size(self):
        """ Estimate of domain size at each vertex as the average of the
            diameters of the circumradius of the tetrahedrons that vertex
            is a part of.

            :rtype: numpy.array
        """
        if self.domain_size is None:
            #coordinates = self.coordinates()
            _ = self.get_vol()

            # Compute the circumradius of the cells
            cr = numpy.zeros((self.tetrahedrons.shape[0]),dtype=float)
            for i in range(len(cr)):
                t_vtx = self.tetrahedrons[i,:]
                # https://en.wikipedia.org/wiki/Tetrahedron#Circumradius
                a = self.distance_between_2_vertices( t_vtx[0], t_vtx[1])
                A = self.distance_between_2_vertices( t_vtx[2], t_vtx[3])
                b = self.distance_between_2_vertices( t_vtx[0], t_vtx[2])
                B = self.distance_between_2_vertices( t_vtx[1], t_vtx[3])
                c = self.distance_between_2_vertices( t_vtx[0], t_vtx[3])
                C = self.distance_between_2_vertices( t_vtx[1], t_vtx[2])
                R = numpy.sqrt( (a*A+b*B+c*C)*(a*A+b*B-c*C)*(a*A-b*B+c*C)*(-a*A+b*B+c*C) ) / (24*self.tetrahedron_vol[i])
                cr[i] = R

            # Compute the mean for each vertex based on all incident cells
            self.domain_size = numpy.zeros((self.vertices.shape[0]),dtype=float)
            count = numpy.zeros((self.vertices.shape[0]),dtype=float)
            for tndx in range(self.tetrahedrons.shape[0]):
                for vndx in self.tetrahedrons[tndx,:]:
                    self.domain_size[vndx] += cr[tndx]
                    count[vndx] += 1
            for vndx in range(len(self.domain_size)):
                self.domain_size[vndx] = self.domain_size[vndx]/count[vndx]

        return self.domain_size

    def distance_between_2_vertices(self, a, b):
        """ Get distance between 2 domain vertices.
            
            :param a: Starting point
            :type a: tuple(float, float, float) or tuple(float, float)

            :param b: Ending point
            :type b: tuple(float, float, float) or tuple(float, float)

            :rtype: float
        """
        return numpy.linalg.norm( self.vertices[a,:]-self.vertices[b,:] )

    def closest_vertex(self, x):
        """ Find the nearest vertex of a given point in the domain.

            :param x: Target source point
            :type x: tuple(float, float, float) or tuple(float, float)

            :rtype: tuple(float, float, float) or tuple(float, float)
        """

        min_dist = None
        min_vtx = None
        for i in range(self.vertices.shape[0]):
            d = numpy.linalg.norm( self.vertices[i,:]-x )
            if d > 0 and (min_dist is None or d < min_dist):
                min_dist = d
                min_vtx = i
        return min_vtx

    def coordinates(self):
        """ Get coordinates within domain.

            :rtype: numpy.array
        """
        return self.vertices

    def get_num_voxels(self):
        """ Get number of voxels in domain.

            :rtype: int
        """

        return self.vertices.shape[0]

    def find_h(self):
        """ Find h value of system.  This value value is based off of \
            the particle which has the greatest distance to \
            its nearest neighbor.

            :rtype: float

        """
        kdtree = KDTree(self.vertices)
        # Detect nearest neighbor distances for all points
        # since each searched point is already included in
        # the tree, we search 2 nearest neighbors, since
        # the first is just the point itself
        distances, indexes = kdtree.query(self.vertices, 2)
        # We only need the distances to the second (non-self) neighbor.
        max_dist = max(distances[:,1])
        h = 2.2*max_dist
        return h

    def get_bounding_box(self):
        """ Get the bounding box of the entire domain.

            :rtype: float | float | float | float | float | float
        """

        xhi=None
        xlo=None
        yhi=None
        ylo=None
        zhi=None
        zlo=None
        for i in range(self.vertices.shape[0]):
            if xhi is None or xhi < self.vertices[i,0]: xhi = self.vertices[i,0]
            if xlo is None or xlo > self.vertices[i,0]: xlo = self.vertices[i,0]
            if yhi is None or yhi < self.vertices[i,1]: yhi = self.vertices[i,1]
            if ylo is None or ylo > self.vertices[i,1]: ylo = self.vertices[i,1]
            if zhi is None or zhi < self.vertices[i,2]: zhi = self.vertices[i,2]
            if zlo is None or zlo > self.vertices[i,2]: zlo = self.vertices[i,2]
        return xhi,xlo,yhi,ylo,zhi,zlo

    def get_vol(self):
        """ Get the total volume of the domain.

            :rtype: float
        """

        if self.vol is None:
               self.calculate_vol()
        return self.vol

    def calculate_vol(self):
        """ Calculate the total volume of the domain.
        """

        self.vol = numpy.zeros((self.vertices.shape[0]),dtype=float)
        self.tetrahedron_vol = numpy.zeros((self.tetrahedrons.shape[0]),dtype=float)
        for t_ndx in range(self.tetrahedrons.shape[0]):
            v1,v2,v3,v4 = self.tetrahedrons[t_ndx]
            a = self.vertices[v1,:]
            b = self.vertices[v2,:]
            c = self.vertices[v3,:]
            d = self.vertices[v4,:]
            #https://en.wikipedia.org/wiki/Tetrahedron#Volume
            t_vol = numpy.abs(numpy.dot( (a-d), numpy.cross( (b-d),(c-d)   ) )/6 )
            self.tetrahedron_vol[t_ndx] = t_vol
            self.vol[v1] += t_vol/4
            self.vol[v2] += t_vol/4
            self.vol[v3] += t_vol/4
            self.vol[v4] += t_vol/4


    def plot_types(self, width=None, height=None, colormap=None, size=5, title=None,
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

        :param return_plotly_figure: Whether or not to return a figure dictionary of data(graph object traces) and layout options
            which may be edited by the user.
        :type return_plotly_figure: bool

        :param use_matplotlib: Whether or not to plot the proprties results using matplotlib.
        :type use_matplotlib: bool
        '''
        from spatialpy.core.Result import _plotly_iterate

        if use_matplotlib:
            width = 6.4 if width in (None, "auto") else width
            height = 4.8 if height in (None, "auto") else height
        else:
            if width in (None, "auto"):
                width = None if width == "auto" else 500
            if height is None:
                height = None if height == "auto" else 500

        if not numpy.count_nonzero(self.vertices[:,1]):
            self.dimensions = 1
        elif not numpy.count_nonzero(self.vertices[:,2]):
            self.dimensions = 2
        else:
            self.dimensions = 3

        if use_matplotlib:
            import matplotlib.pyplot as plt

            if included_types_list is None:
                coords = self.vertices
                type_list = self.type
            else:
                coords = []
                type_list = []
                for i, val in enumerate(self.type):
                    if val in included_types_list:
                        coords.append(self.vertices[i])
                        type_list.append(val)
                coords = numpy.array(coords)

            plt.figure(figsize=(width, height))
            plt.scatter(coords[:,0], coords[:,1], c=type_list, cmap=colormap)
            plt.axis('scaled')
            plt.colorbar()
            if title is not None:
                plt.title(title)
            plt.grid(linestyle='--', linewidth=1)
            plt.plot()
            return

        from plotly.offline import init_notebook_mode, iplot

        types = {}
        for i, val in enumerate(self.type):
            name = "type {}".format(val)

            if included_types_list is None or val in included_types_list:
                if name in types.keys():
                    types[name]['points'].append(self.vertices[i])
                    types[name]['data'].append(self.type[i])
                else:
                    types[name] = {"points":[self.vertices[i]], "data":[self.type[i]]}

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
        else:
            init_notebook_mode(connected=True)
            iplot(fig)


    @classmethod
    def read_xml_mesh(cls, filename):
        """ Read a FEniCS/dolfin style XML mesh file

            :param filename: name of file to read
            :type filename: str

            :rtype: spatialpy.Domain.Domain
        """
        import xml.etree.ElementTree as ET
        root = ET.parse(filename).getroot()
        if not root.tag == 'dolfin': raise DomainError("Not a FEniCS/dolfin xml mesh.")
        mesh = root[0]
        if mesh.tag != 'mesh' or \
           mesh.attrib['celltype'] != 'tetrahedron' or \
           mesh.attrib['dim'] != '3':
            raise DomainError("XML mesh format error")
        #
        vertices = mesh[0]
        cells = mesh[1]
        #vertices
        mesh_vertices = numpy.zeros(( len(vertices), 3), dtype=float)
        for v in vertices:
            mesh_vertices[ int(v.attrib['index']),0] = float(v.attrib['x'])
            mesh_vertices[ int(v.attrib['index']),1] = float(v.attrib['y'])
            mesh_vertices[ int(v.attrib['index']),2] = float(v.attrib['z'])

        # create domain object
        xlim = ( min(mesh_vertices[:,0]) , max(mesh_vertices[:,0]) )
        ylim = ( min(mesh_vertices[:,1]) , max(mesh_vertices[:,1]) )
        zlim = ( min(mesh_vertices[:,2]) , max(mesh_vertices[:,2]) )
        obj = Domain(len(vertices), xlim, ylim, zlim)
        obj.vertices = mesh_vertices

        #tetrahedrons
        obj.tetrahedrons = numpy.zeros(( len(cells), 4), dtype=int)
        for c in cells:
            obj.tetrahedrons[ int(c.attrib['index']),0] = int(c.attrib['v0'])
            obj.tetrahedrons[ int(c.attrib['index']),1] = int(c.attrib['v1'])
            obj.tetrahedrons[ int(c.attrib['index']),2] = int(c.attrib['v2'])
            obj.tetrahedrons[ int(c.attrib['index']),3] = int(c.attrib['v3'])
        # volume
        obj.calculate_vol()
        if not numpy.count_nonzero(obj.vol):
            raise DomainError("Paritcles cannot have 0 volume")
        # set Mass equal to the volume
        obj.mass = obj.vol
        # Calculate density
        obj.rho = obj.mass/obj.vol
        # return model ref
        return obj


    @classmethod
    def import_meshio_object(cls, mesh_obj):
        """ Import a python meshio mesh object.
            
            :param mesh_obj: MeshIO object to import
            :type mesh_obj: meshio.Mesh

            :rtype: spatialpy.Domain.Domain
        """
        # create domain object
        xlim = ( min(mesh_obj.points[:,0]) , max(mesh_obj.points[:,0]) )
        ylim = ( min(mesh_obj.points[:,1]) , max(mesh_obj.points[:,1]) )
        zlim = ( min(mesh_obj.points[:,2]) , max(mesh_obj.points[:,2]) )
        obj = Domain(len(mesh_obj.points), xlim, ylim, zlim)
        #vertices
        obj.vertices = mesh_obj.points
        # triangles
        if 'triangle' in mesh_obj.cells:
            obj.triangles = mesh_obj.cells['triangle']
        #tetrahedrons
        if 'tetra' in mesh_obj.cells:
            obj.tetrahedrons = mesh_obj.cells['tetra']
        # volume
        obj.calculate_vol()
        if not numpy.count_nonzero(obj.vol):
            raise DomainError("Paritcles cannot have 0 volume")
        # set Mass equal to the volume
        obj.mass = obj.vol
        # Calculate density
        obj.rho = obj.mass/obj.vol
        # return model ref
        return obj

    @classmethod
    def read_msh_file(cls, filename):
        """ Read a Gmsh style .msh file
            
            :param filename: Filename of gmsh file
            :type filename: str

            :rtype: spatialpy.Domain.Domain
        """

        try:
            import pygmsh
        except ImportError as e:
            raise DomainError("The python package 'pygmsh' is not installed.")
        try:
            import meshio
        except ImportError as e:
            raise DomainError("The python package 'meshio' is not installed.")

        return cls.import_meshio_object(meshio.msh_io.read(filename))


    def read_stochss_subdomain_file(self, filename):
        """
        Read a .txt file that conains the StochSS v1 spatial subdomain descriptions
        """
        with open(filename,'r') as fd:
            for ln,line in enumerate(fd):
                try:
                    (ndx,type_id) = line.rstrip().split(',')
                    self.type[int(ndx)] = int(type_id)
                except ValueError as e:
                    raise DomainError(f"Could not read in subdomain file, error on line {ln}: {line}")


    @classmethod
    def read_stochss_domain(cls, filename):
        """ Read a StochSS Domain (.domn) file or pull a StochSS Domain from a StochSS Spatial Model (.smdl) file.

            :param filename: name of file to read.
            :type filename: str

            :rtype: spatialpy.Domain.Domain
        """

        try:
            with open(filename, "r") as domain_file:
                domain = json.load(domain_file)
                if "domain" in domain.keys():
                    domain = domain['domain']

            obj = Domain(0, tuple(domain['x_lim']), tuple(domain['y_lim']), tuple(domain['z_lim']),
                        rho0=domain['rho_0'], c0=domain['c_0'], P0=domain['p_0'], gravity=domain['gravity'])

            for particle in domain['particles']:
                # StochSS backward compatability check for rho
                rho = None if "rho" not in particle.keys() else particle['rho']
                # StochSS backward compatability check for c
                c = 0 if "c" not in particle.keys() else particle['c']
                obj.add_point(particle['point'], particle['volume'], particle['mass'],
                              particle['type'], particle['nu'], particle['fixed'], rho=rho, c=c)

            return obj
        except KeyError as e:
            raise DomainError("The file is not a StochSS Domain (.domn) or a StochSS Spatial Model (.smdl).")


    @classmethod
    def create_3D_domain(cls, xlim, ylim, zlim, nx, ny, nz, type_id=1, mass=1.0, nu=1.0, rho=None, c=0, fixed=False, **kwargs):
        """ Create a filled 3D domain

            :param xlim: highest and lowest coordinate in the x dimension
            :type xlim: tuple(float, float)

            :param ylim: highest and lowest coordinate in the y dimension
            :type ylim: tuple(float, float)

            :param zlim: highest and lowest coordinate in the z dimension
            :type zlim: tuple(float, float)

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
            :type rho:

            :param fixed: spatially fixed flag of particles to be created. Defaults to false.
            :type fixed: bool

            :param rho0: background density for the system. Defaults to 1.0
            :type rho0: float

            :param c0: speed of sound for the system. Defaults to 10
            :type c0: float

            :param P0: background pressure for the system. Defaults to 10
            :type P0: float

            :rtype: spatialpy.Domain.Domain
        """
        # Create domain object
        numberparticles = nx*ny*nz
        obj = Domain(numberparticles, xlim, ylim, zlim, **kwargs)
        # Vertices
        obj.vertices = numpy.zeros(( numberparticles, 3), dtype=float)
        x_list = numpy.linspace(xlim[0],xlim[1],nx)
        y_list = numpy.linspace(ylim[0],ylim[1],ny)
        z_list = numpy.linspace(zlim[0],zlim[1],nz)
        ndx = 0
        totalvolume = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0]) * (zlim[1] - zlim[0])
        vol = totalvolume / numberparticles
        if vol < 0:
            raise DomainError("Paritcles cannot have 0 volume")
        for x in x_list:
            for y in y_list:
                for z in z_list:
                    if rho is None:
                        rho = mass / vol
                    obj.vol[ndx] = vol
                    obj.vertices[ndx,0] = x
                    obj.vertices[ndx,1] = y
                    obj.vertices[ndx,2] = z
                    obj.type[ndx] = type_id
                    obj.mass[ndx] = mass
                    obj.nu[ndx] = nu
                    obj.c[ndx] = c
                    obj.rho[ndx] = rho
                    obj.fixed[ndx] = fixed
                    ndx+=1

        # return model ref
        return obj

    @classmethod
    def create_2D_domain(cls, xlim, ylim, nx, ny, type_id=1, mass=1.0, nu=1.0, rho=None, c=0, fixed=False, **kwargs):
        """ Create a filled 2D domain

            :param xlim: highest and lowest coordinate in the x dimension
            :type xlim: tuple(float, float)

            :param ylim: highest and lowest coordinate in the y dimension
            :type ylim: tuple(float, float)

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
            :type rho:

            :param fixed: spatially fixed flag of particles to be created. Defaults to false.
            :type fixed: bool

            :param rho0: background density for the system. Defaults to 1.0
            :type rho0: float

            :param c0: speed of sound for the system. Defaults to 10
            :type c0: float

            :param P0: background pressure for the system. Defaults to 10
            :type P0: float

            :rtype: spatialpy.Domain.Domain
        """
        # Create domain object
        numberparticles = nx*ny
        obj = Domain(numberparticles, xlim, ylim, (0,0), **kwargs)
        # Vertices
        obj.vertices = numpy.zeros(( int(nx)*int(ny), 3), dtype=float)
        x_list = numpy.linspace(xlim[0],xlim[1],nx)
        y_list = numpy.linspace(ylim[0],ylim[1],ny)
        ndx = 0
        totalvolume = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
        vol = totalvolume / numberparticles
        if vol < 0:
            raise DomainError("Paritcles cannot have 0 volume")
        for x in x_list:
            for y in y_list:
                if rho is None:
                    rho = mass / vol
                obj.vol[ndx] = vol
                obj.vertices[ndx,0] = x
                obj.vertices[ndx,1] = y
                obj.vertices[ndx,2] = 0.0
                obj.type[ndx] = type_id
                obj.mass[ndx] = mass
                obj.nu[ndx] = nu
                obj.c[ndx] = c
                obj.rho[ndx] = rho
                obj.fixed[ndx] = fixed
                ndx+=1

        # return model ref
        return obj