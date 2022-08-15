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

#include <time.h>
#include <random>
#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

//#include "binheap.h"
#include "output.h"
#include "particle_system.hpp"
#include "propensities.hpp"
#include "simulate_rdme.hpp"


/**************************************************************************/
namespace Spatialpy{

    void initialize_rdme(ParticleSystem *system, size_t *irN, size_t *jcN,int *prN,
                        size_t *irG,size_t *jcG, unsigned int*u0){

        if(debug_flag){printf("initialize_rdme()\n");fflush(stdout);}
        nsm_core__create(system,irN,jcN,prN,irG,jcG);

        nsm_core__initialize_chem_populations(system, u0);

    }

    /**************************************************************************/
    // This function get called by the main simulation loop.  It advances the
    // state the of RDME system by dt
    void simulate_rdme(ParticleSystem*system,unsigned int step){
        if(!system->static_domain || !system->initialized){
            if(debug_flag) printf("\tnsm_core__initialize_rxn_propensities\n");
            nsm_core__initialize_rxn_propensities(system);
            if(debug_flag) printf("\tnsm_core__initialize_diff_propensities\n");
            nsm_core__initialize_diff_propensities(system);
            if(debug_flag) printf("\tnsm_core__initialize_heap\n");
            if(nsm_core__initialize_heap(system)){
                system->initialized=1;
            }else{
                return; // Don't take step if initialize didn't complete
            }
        }
        if(debug_flag) printf("Simulating RDME for %e seconds\n",system->dt);
        nsm_core__take_step(system, system->dt*step, system->dt);
    }
    /**************************************************************************/
    void destroy_rdme(ParticleSystem*system){
        if(debug_flag) printf("NSM: total # reacton events %lu\n",system->total_reactions);
        if(debug_flag) printf("NSM: total # diffusion events %lu\n",system->total_diffusion);
    }

    void print_current_state(Particle*subvol, ParticleSystem*system){
        unsigned long int i;
        printf("Current state in voxel %i:\n",subvol->id);
        for(i=0;i<system->num_stoch_species;i++){
            printf("xx[%li] = %i\n",i,subvol->xx[i]);
        }
        printf("Neighbors:\n");
        Particle *p2;
        for(auto nn: subvol->neighbors){
            p2 = nn.data;
            printf("%i: nn->D_i_j=%e \n",p2->id,nn.D_i_j);
        }

    }

    /**************************************************************************/
    void nsm_core__create(ParticleSystem*system, size_t *irN, size_t *jcN,int *prN, size_t *irG, size_t *jcG){

        system->irN = irN;
        system->jcN = jcN;
        system->prN = prN;
        system->irG = irG;
        system->jcG = jcG;
        system->total_reactions = 0;
        system->total_diffusion = 0;
        system->initialized = 0;

        Particle*p;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p = &system->particles[i];
            p->srrate = 0;
            p->rrate = (double*) malloc(system->num_stoch_rxns * sizeof(double));
            p->sdrate = 0;
            p->Ddiag = (double*) malloc(system->num_stoch_species * sizeof(double));
        }
    }

    /**************************************************************************/
    void nsm_core__initialize_rxn_propensities(ParticleSystem*system){
        long unsigned int j;
        /* Calculate the propensity for every reaction and every
         subvolume. Store the sum of the reaction intensities in each
         subvolume in srrate. */
        Particle*p;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p = &system->particles[i];
            p->srrate = 0.0;
            for (j = 0; j < system->num_stoch_rxns; j++) {
                double vol = (p->mass / p->rho);
                p->rrate[j] = (*system->stoch_rxn_propensity_functions[j])(p->xx,0.0,vol,p->data_fn,p->type);
                p->srrate += p->rrate[j];
            }
        }
    }

    /**************************************************************************/
    void nsm_core__initialize_diff_propensities(ParticleSystem* system){
        long unsigned int s_ndx;

        Particle *p,*p2;
        double diff_const;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p= &system->particles[i];
            if(p->neighbors.size() == 0){
                p->find_neighbors(system);
            }
            p->sdrate = 0.0;
            for(s_ndx=0; s_ndx<system->num_stoch_species; s_ndx++){
                p->Ddiag[s_ndx] = 0.0;  //Ddiag is sum of (diff_const*n2->D_i_j)
                for(auto n2=p->neighbors.begin(); n2!=p->neighbors.end(); ++n2){
                    p2 = n2->data;
                    diff_const = system->subdomain_diffusion_matrix[s_ndx*system->num_types + (p2->type-1)];
                    p->Ddiag[s_ndx] += diff_const*n2->D_i_j;
                }
                p->sdrate += p->Ddiag[s_ndx] * p->xx[s_ndx];
            }
        }
    }

    /**************************************************************************/
    bool nsm_core__initialize_heap(ParticleSystem*system){
        /** TODO: UPDATE THIS FOR NEW EVENT DATA STRUCTURE **/
        /** TODO: INSERT DR. SANFT'S CODE HERE **/

        /* Calculate times to next event (reaction or diffusion)
         in each subvolume and initialize heap. */

        std::vector<double> propensities(system->particles.size());
        double propensitySum=0.0;
        std::size_t activeChannels=0; // initial number of nonzero propensities (typically improves initial performance)

        long unsigned int p_i = 0 ;
        for(auto p = system->particles.begin(); p!=system->particles.end(); p++){
            propensities[p_i] = p->srrate + p->sdrate;
            propensitySum += propensities[p_i];
            if(propensities[p_i] > 0){
                activeChannels++ ;
            }
            if(p->particle_index != p_i){
                // particles can move around in the sys->particles vector,
                // this ensures that we know where in the vector each particle is
                p->particle_index = p_i;
            }
            p_i++ ;
        }

        double timeOffset = system->current_step * system->dt;
        double endTime  = timeOffset + system->dt;
        if(system->static_domain){
            endTime  = system->nt*system->dt;
        }
        if(system->num_stoch_species > 0){
            if( propensitySum > 0.0 ){
                return system->rdme_event_q.build(propensities, propensitySum,
                                    activeChannels, rng, timeOffset,
                                    endTime );
            }
        }
        return false;

    }

    /**************************************************************************/
    /**************************************************************************/
    void nsm_core__initialize_chem_populations(ParticleSystem*system, unsigned int*u0){
        int num_s = system->num_stoch_species;
        int i = 0 ;
        for(auto p = system->particles.begin(); p!=system->particles.end(); p++){
            p->xx = &u0[i*num_s];
            i++;
        }
    }

    /**************************************************************************/
    /**************************************************************************/
    // Update to use priority queue
    void nsm_core__take_step(ParticleSystem*system, double current_time, double step_size){

        if( system->num_stoch_species == 0 ){ return;}
        double tt = current_time;
        double end_time = current_time + step_size;
        double totrate,cum,rdelta,rrdelta;
        int event,errcode = 0;
        long unsigned int re, spec = 0 ;
        Particle* subvol;
        size_t i,j = 0;
        double rand1,rand2,cum2,old;
        double vol,diff_const;
        Particle *dest_subvol = NULL;

        if(debug_flag>1){
            printf("take_step() tt=%e end_time=%e\n",tt,end_time);
            fflush(stdout);
        }

        std::pair<double,int> timeRxnPair;
        int subvol_index;
        /* Main loop. */
        while(tt <= end_time){
            /* Get the subvolume in which the next event occurred.
             This subvolume is on top of the heap. */

            timeRxnPair = system->rdme_event_q.selectReaction();
            tt = timeRxnPair.first;
            subvol_index = timeRxnPair.second;
            if(subvol_index == -1){ // catch special case of an empty heap
                //printf("ending take_step, subvol_index=%d tt=%e\n",subvol_index,tt);
                //fflush(stdout);
                //debug_flag = 2;
                return;
            }
            subvol = &system->particles[subvol_index];
            vol = (subvol->mass / subvol->rho);
            //if(debug_flag){printf("nsm: tt=%e subvol=%i\n",tt,subvol->id);}
            /* First check if it is a reaction or a diffusion event. */
            totrate = subvol->srrate + subvol->sdrate;


            rand1 = rng() * 1.0 / rng.max();

            if (rand1 <= subvol->srrate/totrate) { // use normalized floating point comparision
                /* Reaction event. */
                event = 0;

                /* a) Determine the reaction re that did occur (direct SSA). */
                double rand_rval = rand1 * subvol->srrate;
                for (re = 0, cum = subvol->rrate[0]; re < system->num_stoch_rxns && rand_rval > cum; re++, cum += subvol->rrate[re]);
                if(re >= system->num_stoch_rxns){
                    if(cum != subvol->srrate){
                        printf("Reaction propensity mismatch in voxel %i. re=%li, srrate[subvol]=%e cum=%e rand_rval=%e\n",subvol->id,re,subvol->srrate,cum,rand_rval);
                        rdelta = 0.0;
                        for (j=0;j<system->num_stoch_rxns; j++) {
                            rdelta += (subvol->rrate[j] = (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type));
                        }
                        subvol->srrate = rdelta;
                    }
                    if(subvol->srrate == 0.0){ continue; }


                    double rand_rval2 = rand1 * subvol->srrate; // sum of propensitiess is not propensity sum, re-roll

                    for (re = 0, cum = subvol->rrate[0]; re < system->num_stoch_rxns && rand_rval2 > cum; re++, cum += subvol->rrate[re]);
                    if(re >= system->num_stoch_rxns){ // failed twice, problems!
                        printf("Propensity sum overflow, rand=%e rand_rval=%e rand_rval2=%e srrate[%i]=%e cum=%e\n",rand1,rand_rval,rand_rval2,subvol->id,subvol->srrate,cum);
                        exit(1);
                    }
                }
                if(debug_flag){printf("nsm: tt=%e subvol=%i type=%i ",tt,subvol->id,subvol->type);}
                if(debug_flag){printf("Rxn %li \n",re);}
                /* b) Update the state of the subvolume subvol and sdrate[subvol]. */
                for (i = system->jcN[re]; i < system->jcN[re+1]; i++) {
                    int prev_val = subvol->xx[system->irN[i]];
                    subvol->xx[system->irN[i]] += system->prN[i];
                    if (subvol->xx[system->irN[i]] < 0){
                        errcode = 1;
                        printf("Negative state detected after reaction %li, subvol %i, species %li at time %e (was %i now %i)\n",re,subvol->id,system->irN[i],tt,prev_val,subvol->xx[system->irN[i]]);

                        print_current_state(subvol,system);
                        exit(1);
                    }
                    subvol->sdrate += subvol->Ddiag[system->irN[i]]*system->prN[i];
                }

                /* c) Recalculate srrate[subvol] using dependency graph. */
                for (i = system->jcG[system->num_stoch_species+re], rdelta = 0.0; i < system->jcG[system->num_stoch_species+re+1]; i++) {
                    old = subvol->rrate[system->irG[i]];
                    j = system->irG[i];
                    rdelta +=
                    (subvol->rrate[j] =
                     (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type)
                     )-old;

                }
                subvol->srrate += rdelta;

                system->total_reactions++; /* counter */
            }
            else {
                /* Diffusion event. */
                event = 1;

                /* a) Determine which species... */
                double diff_rand = rand1 * subvol->sdrate;

                for (spec = 0, cum = subvol->Ddiag[spec]*subvol->xx[spec];
                     spec < system->num_stoch_species && diff_rand > cum;
                     spec++, cum += subvol->Ddiag[spec]*subvol->xx[spec]);
                if(spec >= system->num_stoch_species){
                    // try again, 'cum' is a better estimate of the propensity sum
                    if(cum != subvol->sdrate){
                        printf("Diffusion propensity mismatch in voxel %i. spec=%li, sdrate[subvol]=%e cum=%e diff_rand=%e\n",subvol->id,spec,subvol->sdrate,cum,diff_rand);
                        rdelta = 0.0;
                        for(j = 0; j < system->num_stoch_species; j++){
                            rdelta += subvol->Ddiag[j]*subvol->xx[j];
                        }
                        subvol->sdrate = rdelta;
                    }
                    if(subvol->sdrate == 0.0){ continue; }

                    diff_rand = cum *rand1;
                    for (spec = 0, cum = subvol->Ddiag[spec]*subvol->xx[spec];
                         spec < system->num_stoch_species && diff_rand > cum;
                         spec++, cum += subvol->Ddiag[spec]*subvol->xx[spec]);
                    if(spec >= system->num_stoch_species){
                        spec--;
                        while(subvol->xx[spec] <= 0){
                            spec--;
                            if(spec <=0){
                                printf("Error: diffusion event in voxel %i was selected, but no molecues to move\n",subvol->id);
                                print_current_state(subvol,system);
                                exit(1);
                            }
                        }
                    }
                }


                /* b) and then the direction of diffusion. */
                double r2 = rng() * 1.0 / rng.max();
                rand2 = r2 * subvol->Ddiag[spec];

                /* Search for diffusion direction. */
                Particle*p2;
                cum2 = 0.0;
                for(auto nn=subvol->neighbors.begin(); nn!=subvol->neighbors.end(); ++nn){
                    p2 = nn->data;
                    diff_const = system->subdomain_diffusion_matrix[spec*system->num_types + (p2->type-1)];
                    cum2 += nn->D_i_j * diff_const;
                    if(cum2 > rand2){
                        dest_subvol = p2;
                        break;
                    }
                }
                if(dest_subvol==NULL){ //Diffusion direction overflow
                    // try again, 'cum2' is a better estimate of propensity sum
                    rand2 = r2*cum2;
                    cum2 = 0.0;
                    for(auto nn=subvol->neighbors.begin(); nn!=subvol->neighbors.end(); ++nn){
                        p2 = nn->data;
                        diff_const = system->subdomain_diffusion_matrix[spec*system->num_types + (p2->type-1)];
                        cum2 += nn->D_i_j * diff_const;
                        if(cum2 > rand2){
                            dest_subvol = p2;
                            break;
                        }
                    }
                    if(dest_subvol==NULL){
                        printf("Error: overflow in trying to determine which destination voxel subvol=%i\n",subvol->id);
                        printf("subvol->id=%u ",subvol->id);
                        printf("rand2=%e ",rand2);
                        printf("cum2=%e ",cum2);
                        fflush(stdout);
                        print_current_state(subvol,system);
                        exit(1);
                    }
                }


                /* c) Execute the diffusion event (check for negative elements). */
                subvol->xx[spec]--;
                if (subvol->xx[spec] < 0){
                        errcode = 1;
                        printf("Negative state detected after diffusion, voxel %i -> %i, species %li at time %e\n",subvol->id,dest_subvol->id,spec,tt);
                        print_current_state(subvol,system);
                        exit(1);
                }

                dest_subvol->xx[spec]++;


                if(debug_flag){printf("nsm: tt=%e subvol=%i type=%i ",tt,subvol->id,subvol->type);}
                if(debug_flag){printf("Diff %i->%i",subvol->id,dest_subvol->id);}
                if(debug_flag){
                    printf("    xx[%i]=[",subvol->id);
                    for(i=0;i<system->num_stoch_species;i++){
                        printf("%u ",subvol->xx[i]);
                    }
                    printf("]\n");
                }


                /* Save reaction and diffusion rates. */
                /* Recalculate the reaction rates using dependency graph G. */
                if (system->num_stoch_rxns > 0){
                    for (i = system->jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < system->jcG[spec+1]; i++) {

                        j = system->irG[i];
                        // update this voxel's reactions
                        old = subvol->rrate[j];
                        rdelta +=
                          (subvol->rrate[j] =
                            (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type)
                          )-old;


                        // update dest voxel's reactions
                        old = dest_subvol->rrate[j];
                        rrdelta += (dest_subvol->rrate[j] =
                            (*system->stoch_rxn_propensity_functions[j])(dest_subvol->xx,tt,vol,dest_subvol->data_fn,dest_subvol->type)
                          )-old;
                    }

                    subvol->srrate += rdelta;
                    dest_subvol->srrate += rrdelta;
                }

                /* Adjust diffusion rates. */
                subvol->sdrate -= subvol->Ddiag[spec];
                dest_subvol->sdrate += dest_subvol->Ddiag[spec];

                system->total_diffusion++; /* counter */

            }

            /* Compute time to new event for this subvolume. */
            totrate = subvol->srrate+subvol->sdrate;

            /* Update the heap. */
            system->rdme_event_q.update(subvol_index, totrate, tt, rng);

            /* If it was a diffusion event, also update the other affected
             node. */
            if(event) {
                totrate = dest_subvol->srrate+dest_subvol->sdrate;
                system->rdme_event_q.update(dest_subvol->particle_index, totrate, tt, rng);
            }

            /* Check for error codes. */
            if (errcode) {
                /* Cannot continue. Clear this solution and exit. */
                printf("Exiting due to errcode %i\n",errcode);
                print_current_state(subvol,system);
                exit(1);
            }
        }


    }
}
/**************************************************************************/
