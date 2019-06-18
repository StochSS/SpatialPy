import numpy

class Mesh():
    """ Mesh class for spatial py """


    def __init__(self):
        self.vertices = numpy.zeros((0))
        self.tetrahedrons = numpy.zeros((0,),dtype=int)
        self.vol = None


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
        print("find_h = {0}".format(h))
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
            for t_ndx in range(self.tetrahedrons.shape[0]):
                v1,v2,v3,v4 = self.tetrahedrons[t_ndx]
                a = self.vertices[v1,:]
                b = self.vertices[v2,:]
                c = self.vertices[v3,:]
                d = self.vertices[v4,:]
                #https://en.wikipedia.org/wiki/Tetrahedron#Volume
                t_vol = numpy.abs(numpy.dot( (a-d), numpy.cross( (b-d),(c-d)   ) )/6 )
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
        if not root.tag == 'dolfin': raise MeshException("Not a dolfin mesh.")
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









class MeshException(Exception):
    pass
