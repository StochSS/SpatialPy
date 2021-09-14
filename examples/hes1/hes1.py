try:
    import matplotlib.pyplot as plt
except:
    pass

import os
import numpy as np

import spatialpy

class Nucleus(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return x[0]**2 + x[1]**2 + x[2]**2 <= 3**2
    
class Cytoplasm(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return x[0]**2 + x[1]**2 + x[2]**2 > 3**2

class hes1(spatialpy.Model):
    def __init__(self, model_name="hes1"):
        spatialpy.Model.__init__(self, model_name)

        #Species
        Pf = spatialpy.Species(name="Pf", diffusion_constant=0)
        Po = spatialpy.Species(name="Po", diffusion_constant=0)
        mRNA = spatialpy.Species(name="mRNA", diffusion_constant=6e-1)
        protein = spatialpy.Species(name="protein", diffusion_constant=6e-1)
        self.add_species([Pf, Po, mRNA, protein])
        
        #Restrict to promoter_site
        self.restrict(Pf, 3)
        self.restrict(Po, 3)
        
        #Domain
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.domain = spatialpy.Domain.read_msh_file(basedir + "/mesh/cell.msh")
        
        #Types
        self.set_type(Cytoplasm(), 1)
        self.set_type(Nucleus(), 2)
        # promoter site
        self.domain.type[self.domain.closest_vertex([0]*3)] = 3

        #Parameters
        k1 = spatialpy.Parameter(name="k1", expression=1e9)
        k2 = spatialpy.Parameter(name="k2", expression=0.1)
        alpha_m = spatialpy.Parameter(name="alpha_m", expression=3)
        alpha_m_gamma = spatialpy.Parameter(name="alpha_m_gamma", expression=0.1)
        alpha_p = spatialpy.Parameter(name="alpha_p", expression=1)
        mu_m = spatialpy.Parameter(name="mu_m", expression=0.015)
        mu_p = spatialpy.Parameter(name="mu_p", expression=0.043)
        self.add_parameter([k1, k2, alpha_m, alpha_m_gamma, alpha_p, mu_m, mu_p])


        #Reactions
        R1 = spatialpy.Reaction(name="R1", rate=k1, restrict_to=3,
                                reactants={Pf:1,protein:1}, products={Po:1})
        R2 = spatialpy.Reaction(name="R2", rate=k2, restrict_to=3,
                                reactants={Po:1}, products={Pf:1,protein:1})
        R3 = spatialpy.Reaction(name="R3", rate=alpha_m, restrict_to=3,
                                reactants={Pf:1}, products={Pf:1,mRNA:1})
        R4 = spatialpy.Reaction(name="R4", rate=alpha_m_gamma, restrict_to=3,
                                reactants={Po:1}, products={Po:1,mRNA:1})
        R5 = spatialpy.Reaction(name="R5", rate=alpha_p, restrict_to=1,
                                reactants={mRNA:1}, products={mRNA:1,protein:1})
        R6 = spatialpy.Reaction(name="R6", rate=mu_m, reactants={mRNA:1})
        R7 = spatialpy.Reaction(name="R7", rate=mu_p, reactants={protein:1})
        self.add_reaction([R1, R2, R3, R4, R5, R6, R7])

        #Initail Conditions
        self.add_initial_condition(spatialpy.ScatterInitialCondition(Pf, 1, types=[3]))
        self.add_initial_condition(spatialpy.ScatterInitialCondition(protein, 60, types=[1]))
        self.add_initial_condition(spatialpy.ScatterInitialCondition(mRNA, 10, types=[2]))
        
        self.timespan(range(1200))

if __name__=="__main__":
    model = hes1()
    result = model.run()

    try:
        protein = result.get_species("protein")
        proteinsum = numpy.sum(protein,axis=1)
        plt.plot(model.tspan,proteinsum,'r', label='protein')
        mRNA = result.get_species("mRNA")
        mRNAsum=numpy.sum(mRNA[:],axis=1)
        plt.plot(model.tspan,mRNAsum,'b', label='mRNA')
        plt.legend(loc='best')
        plt.show()
    except:
        pass
