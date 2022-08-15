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

#!/usr/bin/env python3

from spatialpy.Model import SpeciesError
import numpy
import os.path
import spatialpy

class DomainSize(spatialpy.DataFunction):
    def __init__(self, domain):
        spatialpy.DataFunction.__init__(self, name="DomainSize")
        self.domain = domain
        self.h = domain.get_domain_size()
    
    def map(self, point):
        ret = self.h[self.domain.closest_vertex(point)]
        return ret

def create_mincde(model_name="mincde", parameter_values=None):
    """ Model of MinD oscillations in E. Coli, based on the model by Huang. et. al. in """
    model = spatialpy.Model(model_name)

    # Set domain type ids
    model.CYTOSOL = "Cytosol"
    model.MEMBRANE = "Membrane"

    # System constants
    D_m = 1e-14
    D_c = 2.5e-12

    # Make sure that we have the correct path to the mesh file even if we are not executing from the basedir.
    basedir = os.path.dirname(os.path.abspath(__file__))
    domain = spatialpy.Domain.read_xml_mesh(f"{basedir}/data/coli.xml")

    # Set types
    domain.set_properties(spatialpy.GeometryInterior(), model.CYTOSOL)
    domain.set_properties(spatialpy.GeometryExterior(), model.MEMBRATE)

    model.add_domain(domain)

    # Species
    # Restrict to membrane protiens to membrane domain (2)
    MinD_m = spatialpy.Species(name="MinD_m", diffusion_coefficient=D_m, restrict_to=model.MEMBRANE)
    MinD_c_atp = spatialpy.Species(name="MinD_c_atp", diffusion_coefficient=D_c)
    MinD_c_adp = spatialpy.Species(name="MinD_c_adp", diffusion_coefficient=D_c)
    MinD_e = spatialpy.Species(name="MinD_e", diffusion_coefficient=D_c)
    MinDE = spatialpy.Species(name="MinDE", diffusion_coefficient=D_m, restrict_to=model.MEMBRANE)
    model.add_species([MinD_m, MinD_c_atp, MinD_c_adp, MinD_e, MinDE])

    # Parameters
    sigma_d = spatialpy.Parameter(name="sigma_d", expression=2.5e-8)
    sigma_dD = spatialpy.Parameter(name="sigma_dD", expression=0.0016e-18)
    sigma_e = spatialpy.Parameter(name="sigma_e", expression=0.093e-18)
    sigma_de = spatialpy.Parameter(name="sigma_de", expression=0.7)
    sigma_dt = spatialpy.Parameter(name="sigma_dt", expression=1.0)
    model.add_parameter([sigma_d, sigma_dD, sigma_e, sigma_de, sigma_dt])

    # Data function (spatially varying constant)
    model.add_data_function(DomainSize(model.domain))

    # Reactions
    R1 = spatialpy.Reaction(propensity_function="MinD_c_atp * sigma_d / DomainSize", 
                            restrict_to=model.MEMBRANE,
                            reactants={MinD_c_atp: 1}, products={MinD_m: 1})
    R2 = spatialpy.Reaction(rate=sigma_dD,
                            reactants={MinD_c_atp: 1, MinD_m: 1}, products={MinD_m: 2})
    R3 = spatialpy.Reaction(rate=sigma_e,
                            reactants={MinD_m: 1, MinD_e: 1}, products={MinDE: 1})
    R4 = spatialpy.Reaction(rate=sigma_de,
                            reactants={MinDE: 1}, products={MinD_c_adp: 1, MinD_e: 1})
    R5 = spatialpy.Reaction(rate=sigma_dt,
                            reactants={MinD_c_adp: 1}, products={MinD_c_atp: 1})
    R6 = spatialpy.Reaction(rate=sigma_dD,
                            reactants={MinDE: 1, MinD_c_atp: 1}, products={MinD_m: 1, MinDE: 1})
    model.add_reaction([R1, R2, R3, R4, R5, R6])

    # Initial Conditions
    model.add_initial_condition(spatialpy.ScatterInitialCondition(MinD_c_adp, 4500))
    model.add_initial_condition(spatialpy.ScatterInitialCondition(MinD_e, 1575))

    # tspan = spatialpy.TimeSpan(range(200), timestep_size=1)
    tspan = spatialpy.TimeSpan(range(500), timestep_size=1)
    model.timespan(tspan)


if __name__=="__main__":
    """ Dump model to a file. """

    model = create_mincde()
    result = model.run(report_level=1)

    try:
        mindm = result.get_species("MinD_m")
        y_vals = model.mesh.coordinates()[:, 1]
        idx = (y_vals < 1e-6)
        mindmsum = numpy.sum(mindm[:,idx],axis=1)
    except:
        pass
