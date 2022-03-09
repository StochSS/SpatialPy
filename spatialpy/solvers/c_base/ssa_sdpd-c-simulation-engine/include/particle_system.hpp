/**
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2019 - 2022 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef particlesystem_hpp
#define particlesystem_hpp

#include <vector>

#include "ANN/ANN.h" // ANN KD Tree
#include "NRMConstant_v5.hpp"
#include "propensities.hpp"

extern int debug_flag ;

namespace Spatialpy{

    struct Particle;
    struct ParticleSystem;
    struct NeighborNode;
    struct EventNode;

    struct NeighborNode{
        Particle *data;
        double dist;
        double dWdr;
        double D_i_j;
        NeighborNode(Particle *data, double dist, double dWdr, double D_i_j);

        bool operator<(const NeighborNode& n2){
            return dist > n2.dist;
        }
    };

    struct ParticleSystem{
        ParticleSystem(size_t num_types, size_t num_chem_species, size_t num_chem_rxns,
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn);
        ~ParticleSystem();
        int dimension;
        double dt;
        unsigned int nt;
        unsigned int current_step;
        double xlo, xhi, ylo, yhi, zlo, zhi;
        double h;
        double c0;
        double rho0;
        double P0;
        std::vector<Particle> particles;

        // Moved from rdme_t
        NRMConstant_v5 rdme_event_q;
        const size_t *irN;
        const size_t *jcN;
        const int *prN;
        const size_t *irG;
        const size_t *jcG;
        int initialized;
        long int total_reactions;
        long int total_diffusion;

        char boundary_conditions[3];
        bool static_domain;
        size_t num_types;

        ChemRxnFun* chem_rxn_rhs_functions;
        size_t num_chem_species;
        size_t num_chem_rxns;

        PropensityFun* stoch_rxn_propensity_functions;
        size_t num_stoch_species;
        size_t num_stoch_rxns;
        size_t num_data_fn;

        const char * const* species_names;

        const double *subdomain_diffusion_matrix;
        //int *stochic_matrix;
        int *stoichiometric_matrix;
        double* gravity;

        void add_particle(Particle *me);

        ANNkd_tree *kdTree;
        ANNpointArray kdTree_pts;
        bool kdTree_initialized;
    };

/**
system_t* create_system(size_t num_types, size_t num_chem_species, size_t num_chem_rxns,
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn);
particle_t* create_particle(int id);
**/

//bond_t* create_bond(particle p1, particle p2, double k, double rest_distance);

}

#endif //particle_h
