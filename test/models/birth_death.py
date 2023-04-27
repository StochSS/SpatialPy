# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2023 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
''' spatialpy model file for the spatial birth death example. '''
import spatialpy

def create_birth_death():
    ''' Create the spatial birth death model. '''
    model = spatialpy.Model(name='Spatial Birth-Death')

    model.HABITAT = "Habitat"

    domain = spatialpy.Domain.create_2D_domain(
        xlim=(0, 1), ylim=(0, 1), numx=10, numy=10, type_id=model.HABITAT, fixed=True
    )
    model.add_domain(domain)

    rabbits = spatialpy.Species(name='Rabbits', diffusion_coefficient=0.1)
    model.add_species(rabbits)

    init_rabbit_pop = spatialpy.ScatterInitialCondition(species='Rabbits', count=100)
    model.add_initial_condition(init_rabbit_pop)

    k_birth = spatialpy.Parameter(name='k_birth', expression=10)
    k_death = spatialpy.Parameter(name='k_death', expression=0.1)
    model.add_parameter([k_birth, k_death])

    birth = spatialpy.Reaction(
        name='birth', reactants={}, products={"Rabbits":1}, rate="k_birth"
    )
    death = spatialpy.Reaction(
        name='death', reactants={"Rabbits":1}, products={}, rate="k_death"
    )
    model.add_reaction([birth, death])

    tspan = spatialpy.TimeSpan.linspace(t=10, num_points=11, timestep_size=1)
    model.timespan(tspan)
    return model
