/* SSA-SDPD model file, automatically generated by SpatialPy */
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "count_cores.h"
#include "particle.hpp"
#include "propensities.hpp"
#include "simulate.hpp"
#include "simulate_rdme.hpp"
//#include "dSFMT/dSFMT.h"
#include <random>

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
        sys->add_particle(*p);
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


    int main(int argc, char**argv){
        using namespace Spatialpy ;
        //debug_flag = 1;
        //ParticleSystem* system = create_system();
        // Fix particles in space
        //system->static_domain = 1;
        //CONFIG =
        //system->dt = 1;
        //system->nt = 101;
        //system->output_freq = 1;
        //system->h = 0.5;
        //system->rho0 = 1.0;
        //system->c0 = 10;
        //system->P0 = 10;
        // bounding box
        //system->xlo = -5.1;
        //system->xhi = 5.1;
        //system->ylo = -1.1;
        //system->yhi = 1.1;
        //system->zlo = -1.1;
        //system->zhi = 1.1;
        __SYSTEM_CONFIG__
        // create all particles in system
        init_all_particles(system);
        // Setup chemical reaction system
        //initialize_rdme(system, NUM_VOXELS, NUM_SPECIES, NUM_REACTIONS, input_vol, input_sd,
        //                input_data, input_dsize, input_irN, input_jcN, input_prN, input_irG,
        //                input_jcG, input_species_names, input_u0, input_num_subdomain,
        //                input_subdomain_diffusion_matrix);
        __INIT_RDME__

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
        /*
        if(sflag){
            std::mt19937_64 rng(seed) ;
            //dsfmt_init_gen_rand(&dsfmt, seed);
        }else{
            std::mt19937_64 rng((int)time(NULL)+(int)(1e9*clock())) ;
            //dsfmt_init_gen_rand(&dsfmt, (int)time(NULL)+(int)(1e9*clock()));
        }*/

        if(!tflag){
            num_threads = get_num_processors();
            if(num_threads>8){ num_threads=8; }
        }

        run_simulation(num_threads, system);
        exit(0);

}



