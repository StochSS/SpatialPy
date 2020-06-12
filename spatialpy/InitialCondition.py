import numpy
from spatialpy.Model import ModelError


class InitialCondition():
    """ Class used to defined initial conditions in SpatialPy.
        SubClasses must implement the 'apply(model)' method, which 
        direction modifies the model.u0[species,voxel] matrix.
    """

    def apply(self, model):
        raise ModelError("spatialpy.InitialCondition subclasses must implement apply()")

#TODO: implement InitialConditionFromResult()

class PlaceInitialCondition(InitialCondition):
    def __init__(self, species, count, location):
        self.species = species
        self.count = count
        self.location = location

    def apply(self, model):
        spec_name = self.species.name
        for spec_ndx, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species: break
        vtx = model.mesh.closest_vertex(self.location)
        model.u0[spec_ndx, vtx] += self.count

class UniformInitialCondition(InitialCondition):
    def __init__(self, species, count, subdomains=None):
        self.species = species
        self.count = count
        self.subdomains = subdomains

    def apply(self, model):
        spec_name = self.species.name
        for spec_ndx, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species: break

        if self.subdomains is None:
            nvox = model.mesh.get_num_voxels()
            for vtx in range(nvox):
                model.u0[spec_ndx, vtx] += self.count
        else:
            for i in range(model.mesh.get_num_voxels()):
                if model.mesh.sd[i] in self.subdomains:
                    model.u0[spec_ndx, i] += self.count

        

class ScatterInitialCondition(InitialCondition):
    
    def __init__(self, species, count, subdomains=None):
        """ Scatter 'count' of 'species' randomly over the list of subdomains
            (all subdomains if None)."""
        self.species = species
        self.count = count
        self.subdomains = subdomains

    def apply(self, model):
        spec_name = self.species.name
        for spec_ndx, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species: break

        if self.subdomains is None:
            nvox = model.mesh.get_num_voxels()
            for mol in range(self.count):
                vtx = numpy.random.randint(0, nvox)
                model.u0[spec_ndx, vtx] += 1
        else:
            allowed_voxels = []
            for i in range(model.mesh.get_num_voxels()):
                if model.mesh.sd[i] in self.subdomains:
                    allowed_voxels.append(i)
            nvox = len(allowed_voxels)
            if nvox==0: raise ModelError("ScatterInitialCondition has zero voxels to scatter in. Species={0} count={1} subdomains={2}".format(self.species.name, self.count, self.subdomains))
            for mol in range(self.count):
                v_ndx = numpy.random.randint(0, nvox)
                vtx = allowed_voxels[v_ndx]
                model.u0[spec_ndx, vtx] += 1

