/**
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

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

/* SSA-SDPD model file, automatically generated by SpatialPy */
#include <math.h>
#include <memory>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "count_cores.h"
#include "particle_system.hpp"
#include "propensities.hpp"
#include "simulate.hpp"
#include "simulate_rdme.hpp"

int debug_flag ;

namespace Spatialpy{

    /* Species names */
    __DEFINE_SPECIES__

    /* Number of reactions */
    #define NR __NUMBER_OF_REACTIONS__  // Depricated, will remove in future version.
    #define NUM_REACTIONS __NUMBER_OF_REACTIONS__
    #define NUM_SPECIES __NUMBER_OF_SPECIES__
    #define NUM_VOXELS __NUMBER_OF_VOXELS__

    __DATA_FUNCTION_DEFINITIONS__


    /* Parameter definitions */
    __DEFINE_PARAMETERS__

    /* Reaction definitions */
    __DEFINE_REACTIONS__
    /* Deterministic RHS definitions */
    __DEFINE_CHEM_FUNS__

    PropensityFun *ALLOC_propensities(void)
    {
        PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NUM_REACTIONS);
    __DEFINE_PROPFUNS__
        return ptr;
    }

    void FREE_propensities(PropensityFun* ptr)
    {
        free(ptr);
    }

    ChemRxnFun* ALLOC_ChemRxnFun(void){
        ChemRxnFun*ptr = (ChemRxnFun*)malloc(sizeof(ChemRxnFun)*NUM_REACTIONS);
    __DEFINE_CHEM_FUN_INITS__
        return ptr;
    }
    void FREE_ChemRxnFun(ChemRxnFun* ptr){
        free(ptr);
    }


    __INPUT_CONSTANTS__

    //dsfmt_t dsfmt;

    void init_create_particle(ParticleSystem *sys, unsigned int id, double x, double y, double z, int type, double nu, double mass, double rho, int solidTag, int num_chem_species){
        Particle *p  = new Particle(sys, id) ;
        p->x[0] = x;
        p->x[1] = y;
        p->x[2] = z;
        p->id = id;
        p->type = type;
        p->nu = nu;
        p->mass = mass;
        p->rho = rho;
        p->solidTag = solidTag;
        if(num_chem_species > 0){
            for(int i=0;i<num_chem_species;i++){
                p->C[i] = (double) input_u0[id*num_chem_species+i];
            }
        }
        sys->add_particle(p);
    }


    int init_all_particles(ParticleSystem *sys){
        unsigned int id=0;
        __INIT_PARTICLES__
        return id;
    }


    void applyBoundaryConditions(Particle* me, ParticleSystem* system){
    __BOUNDARY_CONDITIONS__
    }

    }/*end namespace*/


    std::mt19937_64 rng ;

    int main(int argc, char**argv){
        using namespace Spatialpy ;

        int num_threads = 1, sflag = 0, tflag = 0, opt;
        long seed;
        while ((opt = getopt(argc, argv, "s:t:")) != -1) {
            switch (opt) {
            case 's':
                seed = atol(optarg);
                sflag = 1;
                break;
            case 't':
                num_threads = atoi(optarg);
                tflag = 1;
                break;
            case '?':
                printf("Usage: %s [OPTION]...\n", argv[0]);
                printf("Example: %s -t 8 -s 1059\n", argv[0]);
                printf("\nOptional arguments:\n");
                printf("  -s Seed value for random number generation.\n");
                printf("  -t Number of threads to use.\n");
                printf("\nIf no arguments are present, seed will be based on the time plus clock and the threads will be set up to 8.\n");
                break;
            }
        }

        rng = sflag ? std::mt19937_64(seed) : std::mt19937_64(std::random_device{}());

        if(!tflag){
            num_threads = get_num_processors();
            if(num_threads>8){ num_threads=8; }
        }
        // Python generated code goes here
        __SYSTEM_CONFIG__
        // create all particles in system
        init_all_particles(system);
        // Setup chemical reaction system
        __INIT_RDME__

        run_simulation(num_threads, system);
        exit(0);
}
