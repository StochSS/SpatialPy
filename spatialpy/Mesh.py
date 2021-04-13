import json
import numpy
from scipy.spatial import KDTree

class Mesh():
    """ Mesh class for spatial py """


    def __init__(self, numpoints, xlim, ylim, zlim, rho0=1.0, c0=10, P0=10, gravity=None):
        self.vertices = numpy.zeros((numpoints, 3), dtype=float)
        self.triangles = None
        self.tetrahedrons = None

        self.on_boundary = None
        self.mesh_size = None
        self.tetrahedron_vol = None

        self.vol = numpy.zeros((numpoints), dtype=float)
        self.mass = numpy.zeros((numpoints), dtype=float)
        self.type = numpy.zeros((numpoints), dtype=int)
        self.nu = numpy.zeros((numpoints), dtype=float)
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
                     f"{pad}Mesh Size: {self.mesh_size}", f"{pad}RHO_0: {self.rho0}", f"{pad}C_0: {self.c0}",
                     f"{pad}P_0: {self.P0}", f"{pad}Gravity: {self.gravity}", f"{pad}X Limit: {self.xlim}",
                     f"{pad}Y Limit: {self.ylim}", f"{pad}Z Limit: {self.zlim}"]
        domain_strs.extend(["", "Paritcles", ""])
        for i, vertex in enumerate(self.vertices):
            v_str = f"{pad}{i+1}: {vertex}\n{pad}   Volume:{self.vol[i]}, Mass: {self.mass[i]}, "
            v_str += f"Type: {self.type[i]}, nu: {self.nu[i]}, Fixed: {self.fixed[i]}"
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

    def add_point(self, x, vol, mass, type, nu, fixed):	
        self.vol = numpy.append(self.vol,vol)	
        self.mass = numpy.append(self.mass,mass)	
        self.type = numpy.append(self.type,type)	
        self.nu = numpy.append(self.nu,nu)	
        self.fixed = numpy.append(self.fixed,fixed)	

        self.vertices = numpy.append(self.vertices, [x], axis=0)

    def find_boundary_points(self):
        if self.on_boundary is None:
            self.on_boundary = numpy.zeros((self.get_num_voxels()), dtype=bool)
            # exterior triangles are part of one-and-only-one tetrahedron
            if self.triangles is None or len(self.triangles) == 0 or len(self.tetrahedrons) == 0:
                return self.on_boundary
            from itertools import combinations
            triangle_in_tetrahedrons_count = {}
            for i in range(self.get_num_voxels()):
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


    def get_mesh_size(self):
        """ Estimate of mesh size at each vertex as the average of the
            diameters of the circumradius of the tetrahedrons that vertex
            is a part of."""
        if self.mesh_size is None:
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
            self.mesh_size = numpy.zeros((self.vertices.shape[0]),dtype=float)
            count = numpy.zeros((self.vertices.shape[0]),dtype=float)
            for tndx in range(self.tetrahedrons.shape[0]):
                for vndx in self.tetrahedrons[tndx,:]:
                    self.mesh_size[vndx] += cr[tndx]
                    count[vndx] += 1
            for vndx in range(len(self.mesh_size)):
                self.mesh_size[vndx] = self.mesh_size[vndx]/count[vndx]

        return self.mesh_size

    def distance_between_2_vertices(self, a, b):
        return numpy.linalg.norm( self.vertices[a,:]-self.vertices[b,:] )

    def closest_vertex(self, x):
        min_dist=None
        min_vtx = None
        for i in range(self.vertices.shape[0]):
            d = numpy.linalg.norm( self.vertices[i,:]-x )
            if d > 0 and (min_dist is None or d < min_dist):
                min_dist = d
                min_vtx = i
        return min_vtx

    def coordinates(self):
        return self.vertices

    def get_num_voxels(self):
        return self.vertices.shape[0]

    def find_h(self):
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
        if self.vol is None:
               self.calculate_vol()
        return self.vol

    def calculate_vol(self):
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



    @classmethod
    def read_xml_mesh(cls, filename):
        """ Read a FEniCS/dolfin style XML mesh file"""
        import xml.etree.ElementTree as ET
        root = ET.parse(filename).getroot()
        if not root.tag == 'dolfin': raise MeshError("Not a FEniCS/dolfin xml mesh.")
        mesh = root[0]
        if mesh.tag != 'mesh' or \
           mesh.attrib['celltype'] != 'tetrahedron' or \
           mesh.attrib['dim'] != '3':
            raise MeshError("XML mesh format error")
        #
        vertices = mesh[0]
        cells = mesh[1]
        #vertices
        mesh_vertices = numpy.zeros(( len(vertices), 3), dtype=float)
        for v in vertices:
            mesh_vertices[ int(v.attrib['index']),0] = float(v.attrib['x'])
            mesh_vertices[ int(v.attrib['index']),1] = float(v.attrib['y'])
            mesh_vertices[ int(v.attrib['index']),2] = float(v.attrib['z'])

        # create mesh object
        xlim = ( min(mesh_vertices[:,0]) , max(mesh_vertices[:,0]) )
        ylim = ( min(mesh_vertices[:,1]) , max(mesh_vertices[:,1]) )
        zlim = ( min(mesh_vertices[:,2]) , max(mesh_vertices[:,2]) )
        obj = Mesh(len(vertices), xlim, ylim, zlim)
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
        # set Mass equal to the volume
        obj.mass = obj.vol
        # return model ref
        return obj


    @classmethod
    def import_meshio_object(cls, mesh_obj):
        """ Import a python meshio mesh object. """
        # create mesh object
        obj = Mesh()
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
        # return model ref
        return obj

    @classmethod
    def read_msh_file(cls, filename):
        """ Read a Gmsh style .msh file """
        try:
            import pygmsh
        except ImportError as e:
            raise MeshError("The python package 'pygmsh' is not installed.")
       # try:
       #     _ = pygmsh.get_gmsh_major_version()
       # except FileNotFoundError as e:
       #     raise MeshError("The command line program 'gmsh' is not installed or is not found in the current PATH")

        try:
            import meshio
        except ImportError as e:
            raise MeshError("The Python package 'meshio' is not installed.")

        return cls.import_meshio_object(meshio.msh_io.read(filename))



    @classmethod
    def read_stochss_domain(cls, filename):
        """
        Read a StochSS Domain (.domn) file or
        pull a StochSS Domain from a StochSS Spatial Model (.smdl) file
        """
        try:
            with open(filename, "r") as domain_file:
                domain = json.load(domain_file)
                if "domain" in domain.keys():
                    domain = domain['domain']

            mesh = Mesh(0, tuple(domain['x_lim']), tuple(domain['y_lim']), tuple(domain['z_lim']),
                        rho0=domain['rho_0'], c0=domain['c_0'], P0=domain['p_0'], gravity=domain['gravity'])

            for particle in domain['particles']:
                mesh.add_point(particle['point'], particle['volume'], particle['mass'],
                               particle['type'], particle['nu'], particle['fixed'])

            return mesh
        except KeyError as e:
            raise MeshError("The file is not a StochSS Domain (.domn) or a StochSS Spatial Model (.smdl).")


    @classmethod
    def create_3D_domain(cls, xlim, ylim, zlim, nx, ny, nz, type_id=1, mass=1.0, nu=1.0, fixed=False, **kwargs):
        """ Create a filled 3D domain
        Args:
            xlim: (tuple) highest and lowest coordinate in the x dimension
            ylim: (tuple) highest and lowest coordinate in the y dimension
            zlim: (tuple) highest and lowest coordinate in the z dimension
            nx: (int) number of particle spacing in the x dimension
            ny: (int) number of particle spacing in the y dimension
            nz: (int) number of particle spacing in the z dimension
            type_id: (int, default: 1) default type ID of particles created to be created
            mass: (float, default: 1.0) default mass of particles created to be created
            nu: (float, default: 1.0) default viscosity of particles created to be created
            fixed: (bool, default: false) spatially fixed flag of particles created to be created
            rho0: (float, default: 1.0) background density for the system
            c0: (float, default: 10) speed of sound for the system
            P0: (float, default: 10) background pressure for the system
        Returns:
            Mesh object
        """
        # Create mesh object
        numberparticles = nx*ny*nz
        obj = Mesh(numberparticles, xlim, ylim, zlim, **kwargs)
        # Vertices
        obj.vertices = numpy.zeros(( numberparticles, 3), dtype=float)
        x_list = numpy.linspace(xlim[0],xlim[1],nx)
        y_list = numpy.linspace(ylim[0],ylim[1],ny)
        z_list = numpy.linspace(zlim[0],zlim[1],nz)
        ndx = 0
        totalvolume = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0]) * (zlim[1] - zlim[0])
        for x in x_list:
            for y in y_list:
                for z in z_list:
                    obj.vol[ndx] = totalvolume / numberparticles
                    obj.vertices[ndx,0] = x
                    obj.vertices[ndx,1] = y
                    obj.vertices[ndx,2] = z
                    obj.type[ndx] = type_id
                    obj.mass[ndx] = mass
                    obj.nu[ndx] = nu
                    obj.fixed[ndx] = fixed
                    ndx+=1

        # return model ref
        return obj

    @classmethod
    def create_2D_domain(cls, xlim, ylim, nx, ny, type_id=1, mass=1.0, nu=1.0, fixed=False, **kwargs):
        """ Create a filled 2D domain
        Args:
            xlim: (tuple) highest and lowest coordinate in the x dimension
            ylim: (tuple) highest and lowest coordinate in the y dimension
            nx: (int) number of particle spacing in the x dimension
            ny: (int) number of particle spacing in the y dimension
            type_id: (int, default: 1) default type ID of particles created to be created
            mass: (float, default: 1.0) default mass of particles created to be created
            nu: (float, default: 1.0) default viscosity of particles created to be created
            fixed: (bool, default: false) spatially fixed flag of particles created to be created
            rho0: (float, default: 1.0) background density for the system
            c0: (float, default: 10) speed of sound for the system
            P0: (float, default: 10) background pressure for the system
        Returns:
            Mesh object
        """
        # Create mesh object
        numberparticles = nx*ny
        obj = Mesh(numberparticles, xlim, ylim, (0,0), **kwargs)
        # Vertices
        obj.vertices = numpy.zeros(( int(nx)*int(ny), 3), dtype=float)
        x_list = numpy.linspace(xlim[0],xlim[1],nx)
        y_list = numpy.linspace(ylim[0],ylim[1],ny)
        ndx = 0
        totalvolume = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
        #print("totalvolume",totalvolume)
        for x in x_list:
            for y in y_list:
                obj.vol[ndx] = totalvolume / numberparticles
                obj.vertices[ndx,0] = x
                obj.vertices[ndx,1] = y
                obj.vertices[ndx,2] = 0.0
                obj.type[ndx] = type_id
                obj.mass[ndx] = mass
                obj.nu[ndx] = nu
                obj.fixed[ndx] = fixed
                ndx+=1

        # return model ref
        return obj



class MeshError(Exception):
    pass
