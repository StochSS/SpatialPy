import numpy

class Mesh():
    """ Mesh class for spatial py """


    def __init__(self):
        self.vertices = numpy.zeros((0))
        self.tetrahedrons = numpy.zeros((0,),dtype=int)


    def coordinates(self):
        return self.vertices

    def get_num_voxels(self):
        return self.vertices.shape[0]


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
        obj.tetrahedrons = numpy.zeros(( len(cells), 3), dtype=int)
        for c in cells:
            obj.tetrahedrons[ int(c.attrib['index']),0] = int(c.attrib['v0'])
            obj.tetrahedrons[ int(c.attrib['index']),1] = int(c.attrib['v1'])
            obj.tetrahedrons[ int(c.attrib['index']),2] = int(c.attrib['v2'])
            obj.tetrahedrons[ int(c.attrib['index']),2] = int(c.attrib['v3'])
        # return model ref
        return obj









class MeshException(Exception):
    pass
