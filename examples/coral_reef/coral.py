#!/usr/bin/env python
import os
import numpy
import sys
from collections import OrderedDict
sys.path.append("../..")
import spatialpy

class CoralReef(spatialpy.Model):

    def __init__(self, name="coral_reef", D_c=1.0, D_m=1.0, version=1):
        spatialpy.Model.__init__(self, name)

        # Species
        Coral = spatialpy.Species(name="Coral",diffusion_constant=0.0)
        Coral_m = spatialpy.Species(name="Coral_m",diffusion_constant=D_c)
        MA = spatialpy.Species(name="MA", diffusion_constant=0.0)
        MA_m = spatialpy.Species(name="MA_m", diffusion_constant=D_m)
        Turf = spatialpy.Species(name="Turf", diffusion_constant=0.0)
        self.add_species([Coral, MA, Coral_m, MA_m, Turf])

        # Parameters
        phi_c = spatialpy.Parameter(name="phi_c", expression=0.0011) #1/year
        phi_m = spatialpy.Parameter(name="phi_m", expression=0.001) #1/year
        g_tc = spatialpy.Parameter(name="g_tc", expression=0.1) #1/year
        g_tm = spatialpy.Parameter(name="g_tm", expression=0.2) #1/year
        Gamma = spatialpy.Parameter(name="Gamma", expression=0.05)
        dc = spatialpy.Parameter(name="dc", expression=0.05) #1/year
        dm = spatialpy.Parameter(name="dm", expression=1.0) #1/year
        #dm = spatialpy.Parameter(name="dm", expression=0.2) #1/year
        phi_g = spatialpy.Parameter(name="psi_g", expression=0.0)
        # Death rate of mobile propgules.  Combine with diffusion to determine spread.
        mu_c = spatialpy.Parameter(name="mu_c", expression=1.0) #1/year
        mu_m = spatialpy.Parameter(name="mu_m", expression=1.0) #1/year
        # mobile propogules destroyed by estabilished 
        alpha_c = spatialpy.Parameter(name="alpha_c", expression=0.1) #1/year
        alpha_m = spatialpy.Parameter(name="alpha_m", expression=0.5) #1/year
        # Production of mobile propogules
        R_c = spatialpy.Parameter(name="R_c", expression=1.0) #1/year
        R_m = spatialpy.Parameter(name="R_m", expression=1.0) #1/year

        self.add_parameter([phi_c, phi_m, g_tc, g_tm, Gamma, dc, dm, phi_g, mu_c, mu_m, alpha_c, alpha_m, R_c, R_m])
                

        # Reactions:
        # C -> T : dc
        self.add_reaction(spatialpy.Reaction(name="R3", reactants={Coral:1}, products={Turf:1}, rate=dc))
        # MA -> T : dm
        self.add_reaction(spatialpy.Reaction(name="R4", reactants={MA:1}, products={Turf:1}, rate=dm))
        # T + C_m -> C : phi_c
        self.add_reaction(spatialpy.Reaction(name="R5", reactants={Turf:1, Coral_m:1}, products={Coral:1}, rate=phi_c))
        # T + MA_m -> MA : phi_m
        self.add_reaction(spatialpy.Reaction(name="R6", reactants={Turf:1, MA_m:1}, products={MA:1}, rate=phi_m))
        # C + T -> 2C : g_tc * exp(-1.0 * psi_g * MA / 100)
        self.add_reaction(spatialpy.Reaction(name="R7", reactants={Turf:1, Coral:1}, products={Coral:2}, propensity_function="g_tc*Turf*Coral*exp(-1.0 * psi_g * MA / Space_per_voxel)/vol"))
        # MA + T -> 2MA : g_tm
        self.add_reaction(spatialpy.Reaction(name="R8", reactants={Turf:1, MA:1}, products={MA:2}, rate=g_tm))
        # C + MA -> 2MA : Gamma * g_tm
        self.add_reaction(spatialpy.Reaction(name="R9", reactants={Coral:1, MA:1}, products={MA:2}, propensity_function="g_tm*Gamma*Coral*MA/vol"))
        # C -> C + C_m : R_c
        self.add_reaction(spatialpy.Reaction(name="R10", reactants={Coral:1}, products={Coral:1, Coral_m:1}, rate=R_c))
        # MA -> MA + MA_m : R_m
        self.add_reaction(spatialpy.Reaction(name="R11", reactants={MA:1}, products={MA:1, MA_m:1}, rate=R_m))
        # C_m -> 0 : mu_c
        self.add_reaction(spatialpy.Reaction(name="R12", reactants={Coral_m:1}, products={}, rate=mu_c))
        # MA_m -> 0 : mu_m
        self.add_reaction(spatialpy.Reaction(name="R13", reactants={MA_m:1},  products={}, rate=mu_m))
        # MA + C_m -> MA : alpha_c
        self.add_reaction(spatialpy.Reaction(name="R14", reactants={MA:1, Coral_m:1},  products={MA:1}, rate=alpha_c))
        # C + MA_m -> C : alpha_m
        self.add_reaction(spatialpy.Reaction(name="R15", reactants={Coral:1, MA_m:1},  products={Coral:1}, rate=alpha_m))
 
        

        # A unit square
        # each grid point is 10cm x 10cm, domain is 5m x 5m
        self.mesh = spatialpy.MODELMesh.generate_square_mesh(L=5, nx=50, ny=50, periodic=True)

        Space_per_voxel = 10
        self.add_parameter(spatialpy.Parameter(name="Space_per_voxel", expression=Space_per_voxel)) #1/year
        
                          
        if True:
            # Start with two colonys
            self.set_initial_condition_distribute_uniformly({Turf:Space_per_voxel})
            self.set_initial_condition_place_near({Coral:Space_per_voxel}, point=[1,1])
            self.set_initial_condition_place_near({Turf:0}, point=[1,1])
            self.set_initial_condition_place_near({MA:Space_per_voxel}, point=[4,4])
            self.set_initial_condition_place_near({Turf:0}, point=[4,4])
        else:
            # Every voxel is the same
            self.set_initial_condition_distribute_uniformly({Turf:0})
            self.set_initial_condition_distribute_uniformly({Coral:Space_per_voxel-1})
            self.set_initial_condition_distribute_uniformly({MA:1})


        for vndx in range(self.u0.shape[1]):
            tot = 0
            for sndx, sname in enumerate(self.listOfSpecies):
                tot += self.u0[sndx][vndx]
            if tot > 100:
                for sndx, sname in enumerate(self.listOfSpecies):
                    print("u0[{0}][{1}] = {2}").format(sname, vndx, self.u0[sndx][vndx])

        #self.timespan(numpy.linspace(0,500,501)) #500 years
        #self.timespan(numpy.linspace(0,5,72)) #5 years, by months
        self.timespan(numpy.linspace(0,11,66)) #10 years, by 2 months


if __name__ == "__main__":
    model = CoralReef()
    result = model.run(report_level=1)

    x_vals = model.mesh.coordinates()[:, 0]
    y_vals = model.mesh.coordinates()[:, 1]
    C_vals = result.get_species("Coral")
    MA_vals = result.get_species("MA")
    Turf_vals = result.get_species("Turf")
    num_vox = len(x_vals)    

    print(C_vals)
    print(MA_vals)
    print(Turf_vals)
  
