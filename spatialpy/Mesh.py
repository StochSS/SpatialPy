import numpy

class Mesh():
    """ Mesh class for spatial py """


    def __init__(self):
        self.vertices = numpy.zeros((0))
        self.triangles = numpy.zeros((0),dtype=int)
        self.tetrahedrons = numpy.zeros((0,),dtype=int)
        self.on_boundary = None
        self.vol = None
        self.mesh_size = None
        self.tetrahedron_vol = None

    def find_boundary_points(self):
        if self.on_boundary is None:
            self.on_boundary = numpy.zeros((self.get_num_voxels()),dtype=bool)
            # exterior triangles are part of one-and-only-one tetrahedron
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
        max_dist = None
        #print("find_h")
        for i in range(self.vertices.shape[0]):
            d = self.dist_to_closest_neighbor(i)
            #print("\tdist_to_closest_neighbor({0})={1}".format(i,d))
            if max_dist is None or d > max_dist:
                max_dist = d
        h = 2.2*max_dist
        #print("find_h = {0}".format(h))
        return h


    def dist_to_closest_neighbor(self, v_ndx):
        min_dist=None
        for i in range(self.vertices.shape[0]):
            if i==v_ndx: continue
            d = numpy.linalg.norm( self.vertices[i,:]-self.vertices[v_ndx,:] )
            if d > 0 and (min_dist is None or d < min_dist):
                min_dist = d
        return min_dist

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
                
        return self.vol            


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
            raise MeshException("XML mesh format error")
        #
        vertices = mesh[0]
        cells = mesh[1]
        # create mesh object
        obj = Mesh()
        #vertices
        obj.vertices = numpy.zeros(( len(vertices), 3), dtype=float)
        for v in vertices:
            obj.vertices[ int(v.attrib['index']),0] = float(v.attrib['x'])
            obj.vertices[ int(v.attrib['index']),1] = float(v.attrib['y'])
            obj.vertices[ int(v.attrib['index']),2] = float(v.attrib['z'])
        #tetrahedrons
        obj.tetrahedrons = numpy.zeros(( len(cells), 4), dtype=int)
        for c in cells:
            obj.tetrahedrons[ int(c.attrib['index']),0] = int(c.attrib['v0'])
            obj.tetrahedrons[ int(c.attrib['index']),1] = int(c.attrib['v1'])
            obj.tetrahedrons[ int(c.attrib['index']),2] = int(c.attrib['v2'])
            obj.tetrahedrons[ int(c.attrib['index']),3] = int(c.attrib['v3'])
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
        # return model ref
        return obj

    @classmethod
    def read_msh_file(cls, filename):
        """ Read a Gmsh style .msh file """
        try:
            import pygmsh
        except ImportError as e:
            raise MeshError("The python package 'pygmsh' is not installed.")
        try:
            _ = pygmsh.get_gmsh_major_version()
        except FileNotFoundError as e:
            raise MeshError("The command line program 'gmsh' is not installed or is not found in the current PATH")

        try:
            import meshio
        except ImportError as e:
            raise MeshError("The python package 'meshio' is not istaled.")
        
        return cls.import_meshio_object(meshio.msh_io.read(filename))







class MeshError(Exception):
    pass
