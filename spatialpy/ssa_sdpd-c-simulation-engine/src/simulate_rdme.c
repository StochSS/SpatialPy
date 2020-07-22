#include "linked_list.h"
#include <time.h>
#include "output.h"
#include "particle.h"
#include "simulate_rdme.h"
#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "propensities.h"
#include "binheap.h"


/**************************************************************************/
void initialize_rdme(system_t*system, const int Ncells, const int Mspecies, const int Mreactions, 
                        size_t *irN, size_t *jcN,int *prN,size_t *irG,size_t *jcG,
                        const char* const species_names[], const unsigned int*u0,
                        const int num_subdomains, const double*subdomain_diffusion_matrix
                        ){
    if(debug_flag){printf("*************** initialize_rdme ******************\n");}
    nsm_core__create(system,irN,jcN,prN,irG,jcG);

    nsm_core__initialize_chem_populations(system, u0);

}

/**************************************************************************/
// This function get called by the main simulation loop.  It advances the
// state the of RDME system by dt
void simulate_rdme(system_t*system,unsigned int step){
    rdme_t*rdme = system->rdme;
    if(rdme == NULL){
        return;
    }
    if(!system->static_domain || !rdme->initialized){
// All of below is replaced by code in find_neighbor()
//            // if the  domain is not static, rebuild the diffusion matrix after movement
//        if(!rdme->initialized){
//            if(debug_flag) printf("Building diffusion matrix\n");
//            if(debug_flag) printf("\tnsm_core__build_diffusion_matrix\n");
//            nsm_core__build_diffusion_matrix(rdme,system);
            rdme->initialized=1;
//        }else{
//            if(debug_flag) printf("Rebuilding diffusion matrix\n");
//            if(debug_flag) printf("\tnsm_core__destroy_diffusion_matrix\n");
//            nsm_core__destroy_diffusion_matrix(rdme);
//            if(debug_flag) printf("\tnsm_core__build_diffusion_matrix\n");
//            nsm_core__build_diffusion_matrix(rdme,system);
//        }
        if(debug_flag) printf("\tnsm_core__initialize_rxn_propensities\n");
        nsm_core__initialize_rxn_propensities(rdme);
        if(debug_flag) printf("\tnsm_core__initialize_diff_propensities\n");
        nsm_core__initialize_diff_propensities(rdme);
        if(debug_flag) printf("\tnsm_core__initialize_heap\n");
        nsm_core__initialize_heap(system);
    }
    if(debug_flag) printf("Simulating RDME for %e seconds\n",system->dt);
    nsm_core__take_step(system, system->dt*step, system->dt);
}
/**************************************************************************/
void destroy_rdme(system_t*system){
    if(system->rdme == NULL){
        return;
    }
    if(debug_flag) printf("NSM: total # reacton events %lu\n",system->rdme->total_reactions);
    if(debug_flag) printf("NSM: total # diffusion events %lu\n",system->rdme->total_diffusion);
    nsm_core__destroy(system->rdme);
}




//===================================================
// Adapted from PyURDME's nsmcore.c
//===================================================

void print_current_state(int subvol, unsigned int*xx,const size_t Mspecies){
    int i;
    printf("Current state in voxel %i:\n",subvol);
    for(i=0;i<Mspecies;i++){
        printf("xx[%i] = %i\n",i,xx[subvol*Mspecies+i]);
    }
}

/*void nsm_core(const size_t *irD,const size_t *jcD,const double *prD,
              const size_t *irN,const size_t *jcN,const int *prN,
              const size_t *irG,const size_t *jcG,
          const double *vol,
          const int *sd,
          const double *data,
              const size_t Ncells,
              const size_t Mspecies,
          const size_t Mreactions,
          const size_t dsize,
          int report_level
              int *xx,
          double end_time)
*/
/* Specification of the inputs:

 Ncells
 Number of subvolumes.

 Mspecies
 Number of species.

 Hence Ndofs = Ncells*Mspecies.

 Mreactions
 Total number of reactions.

 dsize
 Size of data vector sent to propensities.

 end_time
 length of the simulation.  This function simulates the RDME on the
 time interval [0, end_time].

 report_level
 The desired degree of feedback during simulations. 0, 1, and 2 are
 currently supported options.

 Diffusion matrix D. Double sparse (Ndofs X Ndofs).
 Macroscopic diffusion matrix. D(i,j) is the diffusion rate from dof #j to
 dof #i. This matrix uses the CSR-format and not CSC because fast access to
 rows is needed.

 State vector 'xx'. Integer (Mspecies X Ncells).
 Gives the initial copy number of the species in each subvolume.
 !!! This is also the output.

 Stochiometric matrix N. Integer sparse (Mspecies X Nreactions).
 N(:,j) describes how reaction j changes the number of species.

 Dependency graph G. Integer sparse (Mreactions X Mspecies+Mreactions).
 G(i,Mspecies+j) is non-zero if executing reaction j means that reaction i
 needs to be re-evaluated. The first Mspecies columns of G similarily cover
 diffusion events.
 vol. Double vector (length Ncells).
 vol[i] gives the volume of cell #i.

 data. Double matrix (dsize X Ncells).
 Generalized data matrix, data(:,j) gives a data vector for cell #j.

 sd. Integer vector (length Ncells).
 Subdomain number. sd[i] is the subdomain of cell #i. The vector sd can also
 be used to separate boundaries, line segments and points.

 Format of sparse matrices:
 G, N and S are sparse matrices in compressed column format (CCS). D is sparse
 but in compressed row format (CRS), or equivalently, a transposed matrix in
 CCS format.
 jcD, irD, prD (double *)
 jcN, irN, prN (int *)
 jcG, irG (int *)

 Propensities:
 a vector of function pointers (length Mreactions) is input by
 linking with the prototypes in propensities.h and function
 definitions in a user-specified .c-file. The type of this vector is
 PropensityFun which defines the input to a property function. See
 propensities.h for more details.

 Ordering of the dofs:
 Dof #i is located in cell #(i/Mspecies), and the dofs located in
 cell #j is u0(:,j). Thus, u0 is understood as a matrix of size
 Mspecies X Ncells.
 */

/**************************************************************************/
void nsm_core__create(system_t*system, size_t *irN, size_t *jcN,int *prN, size_t *irG, size_t *jcG){
    /* Create the RDME object */
    rdme_t* rdme = (rdme_t*) malloc(sizeof(rdme_t));

    rdme->irN = irN;
    rdme->jcN = jcN;
    rdme->prN = prN;
    rdme->irG = irG;
    rdme->jcG = jcG;
    rdme->total_reactions = 0;
    rdme->total_diffusion = 0;
    rdme->initialized = 0;

    rdme->heap = create_ordered_list();



    node*n;
    particle_t*p;
    int i=0;
    for(n=system->particle_list->head; n!=NULL; n=n->next){
        p = n->data;
        p->rdme = (rdme_voxel_t*) malloc(sizeof(rdme_voxel_t));
        p->rdme->srrate = 0;
        p->rdme->rrate = (double*) malloc(system->num_stoch_species * sizeof(double));
        p->rdme->sdrate = 0;
        p->rdme->Ddiag = (double*) malloc(system->num_stoch_species * sizeof(double));

        ordered_list_add( rdme->heap, p );

    }



    /* Create reaction rate matrix (Mreactions X Ncells) and total rate
     vector. In rrate we store all propensities for chemical rections,
     and in srrate the sum of propensities in every subvolume. */
    //rdme->rrate = (double *)malloc(system->num_stoch_rxn*rdme->Ncells*sizeof(double));
    //rdme->srrate = (double *)malloc(rdme->Ncells*sizeof(double));

    //nsm_core__initialize_rxn_propensities(rdme);

    /* Total diffusion rate vector (length Mcells). It will hold
     the total diffusion rates in each subvolume. */
    //rdme->sdrate = (double *)malloc(rdme->Ncells*sizeof(double));
    /* The diagonal value of the D-matrix is used frequently. For
     efficiency, we store the negative of D's diagonal in Ddiag. */
    //rdme->Ddiag = (double *)malloc(rdme->Ndofs*sizeof(double));

    //nsm_core__initialize_diff_propensities(rdme);

    /* Create binary (min)heap. */
    //rdme->rtimes = (double *)malloc(rdme->Ncells*sizeof(double));
    //rdme->node = (int *)malloc(rdme->Ncells*sizeof(int));
    //rdme->heap = (int *)malloc(rdme->Ncells*sizeof(int));


    /* return rdme structure */
   //return rdme;
   system->rdme = rdme;
}

/**************************************************************************/
void nsm_core__initialize_rxn_propensities(rdme_t* rdme){
    int i,j;
    /* Calculate the propensity for every reaction and every
     subvolume. Store the sum of the reaction intensities in each
     subvolume in srrate. */
    node*n;
    particle_t*p;
    //int i=0;
    for(n=system->particle_list->head; n!=NULL; n=n->next){
    //    p = n->data;
    ///    rdme->sd[p->id] = p->type;
    //    rdme->vol[p->id] = p->mass / p->rho;
    //}
    //for (i = 0; i < rdme->Ncells; i++) {
        p->rdme->srrate = 0.0;
        for (j = 0; j < system->num_stoch_species; j++) {
            //rrate[i*Mreactions+j] =
            //(*rfun[j])(&xx[i*Mspecies],tt,vol[i],&data[i*dsize],sd[i],i,xx,irK,jcK,prK);
            //srrate[i] += rrate[i*Mreactions+j];
            //rdme->rrate[i*system->num_stoch_rxn+j] = (*rdme->rfun[j])(&rdme->xx[i*system->num_stoch_species],0.0,rdme->vol[i],&rdme->data[i*rdme->dsize],rdme->sd[i]);
            //rdme->srrate[i] += rdme->rrate[i*system->num_stoch_rxn+j];
            double vol = (p->mass / p->rho);
            p->rdme->rrate[j] = (*system->stoch_rxn_propensity_functions[j])(p->xx,0.0,vol,p->data_fn,p->type);
            p->rdme->srrate += rdme->rrate[j];
        }
    }
}

/**************************************************************************/
void nsm_core__initialize_diff_propensities(system_t* system){
    int spec;

//    for (i = 0; i < rdme->Ndofs; i++) {
//        rdme->Ddiag[i] = 0.0;
//        for (j = rdme->jcD[i]; j < rdme->jcD[i+1]; j++)
//        if (rdme->irD[j] == i) rdme->Ddiag[i] = -1*rdme->prD[j];
//    }
//
//    /* Calculate the total diffusion rate for each subvolume. */
//    for(i = 0; i < rdme->Ncells; i++) {
//        rdme->sdrate[i] = 0.0;
//        for(j = 0; j < system->num_stoch_species; j++)
//        rdme->sdrate[i] += rdme->Ddiag[i*system->num_stoch_species+j]*rdme->xx[i*system->num_stoch_species+j];
//    }
    node*n;
    neighbor_node_t*n2;
    particle_t*p;
    for(n=system->particle_list->head; n!=NULL; n=n->next){
        p=n->data;
        if(p->neighbors->count == 0){
            find_neighbors(p, system);
        }
        p->rdme->sdrate = 0.0;
        for(s_ndx=0; s_ndx<system->num_stoch_species; s_ndx++){
            p->rdme->Ddiag[s_ndx] = 0.0;  //Ddiag is sum of (diff_const*n2->D_i_j)
            for(n2=p->neighbors->head; n2!=NULL; n2=n2->next){
                p2 = n2->data;
                diff_const = system->subdomain_diffusion_matrix[spec*system->num_subdomains + (p2->type-1)];
                p->rdme->Ddiag[s_ndx] += diff_const*n2->D_i_j;
            }
            p->rdme->sdrate += p->rdme->Ddiag[s_ndx] * p->xx[s_ndx];
        }
    }

}

/**************************************************************************/
void nsm_core__initialize_heap(system_t*system){
    int i;
    rdme_t* rdme = system->rdme;
    /* Calculate times to next event (reaction or diffusion)
     in each subvolume and initialize heap. */
    node_ordered_node_t*n;
    particle_t*p;
    for(n=rdme->heap->head; n!=NULL; n=n->next){
    //for (i = 0; i < rdme->Ncells; i++) {
        //rdme->rtimes[i] = -log(1.0-drand48())/(rdme->srrate[i]+rdme->sdrate[i]);
        //rdme->heap[i] = rdme->node[i] = i;
        p = n->data;
        n->tt = -log(1.0-drand48())/(p->rdme->srrate+p->rdme->sdrate);
    }
    //initialize_heap(rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);
    ordered_list_sort(rdme->heap);
}

/**************************************************************************/
void nsm_core__destroy(rdme_t*rdme){
    FREE_propensities(rdme->rfun);
    free(rdme->xx);
    free(rdme->heap);
    free(rdme->node);
    free(rdme->rtimes);
    free(rdme->Ddiag);
    free(rdme->sdrate);
    free(rdme->srrate);
    free(rdme->rrate);
    free(rdme);
}

/**************************************************************************/
void nsm_core__initialize_chem_populations(rdme_t* rdme, const unsigned int*u0){
    /* Set xx to the initial state. xx will always hold the current solution. */
    //printf("malloc Ndofs = %li\n",rdme->Ndofs);
    //TODO
    rdme->xx = (unsigned int *)malloc(rdme->Ndofs*sizeof(unsigned int));
    memcpy(rdme->xx,u0,rdme->Ndofs*sizeof(unsigned int));
    //printf("       Ndofs = %li\n",rdme->Ndofs);
    //printf("xx = [ ");
    //int i;
    //for(i=0;i<rdme->Ndofs;i++){
    //    printf("%u ",rdme->xx[i]);
    //}
    //printf("]\n");
}


/**************************************************************************/
//void nsm_core__build_diffusion_matrix(rdme_t*rdme,system_t*system){
//    if(debug_flag){ printf("*************** build_diffusion_matrix ***************\n");fflush(stdout);}
//    double off_diag_sum,diff_const,dist2;
//    node *n;
//    neighbor_node_t*n2;
//    particle_t *p1,*p2;
//    int s_ndx;
//    double D_i_j;
//    double ih,ihsq,wfd;
//    double h = system->h;
//
    //size_t jcD_length = rdme->Ncells + 1;
    //size_t irD_length = 0;
    //size_t prD_length = 0;
    // find total length of jc & pr arrays: O(n)
    //for(n=system->particle_list->head; n!=NULL; n=n->next){
    //    p1 = n->data;
    //    if(p1->neighbors->count == 0){
    //        //if(debug_flag){printf("find_neighbors(%i)\n",p1->id);}
    //        find_neighbors(p1, system);
    //    }
    //    //if(debug_flag){printf("node %i # neighbors %i\n",p1->id,p1->neighbors->count);}
    //    //irD_length += (p1->neighbors->count + 1);
    //}
    //prD_length = irD_length;
    //if(debug_flag){printf("irD_length= %li\n",irD_length);fflush(stdout);}
    //if(debug_flag){printf("jcD_length= %li\n",jcD_length);fflush(stdout);}
    // allocate space for each array
    //printf("MALLOC rdme->irD [%li]\n",irD_length*system->num_stoch_species);
    //rdme->irD = (size_t*) malloc(sizeof(size_t)*irD_length*system->num_stoch_species);
    //size_t irD_ndx = 0;
    //printf("MALLOC rdme->jcD [%li]\n",jcD_length*system->num_stoch_species);
    //rdme->jcD = (size_t*) malloc(sizeof(size_t)*jcD_length*system->num_stoch_species);
    //size_t jcD_ndx = 0;
    //rdme->jcD[jcD_ndx++] = 0;
    //printf("MALLOC rdme->prD [%li]\n",prD_length*system->num_stoch_species);
    //rdme->prD = (double*) malloc(sizeof(double)*prD_length*system->num_stoch_species);
    //size_t prD_ndx = 0;
    // for each particle p, look at each neighbor p2
//    for(n=system->particle_list->head; n!=NULL; n=n->next){
//        p1 = n->data;
//        if(p1->neighbors->count == 0){
//            find_neighbors(p1, system);
//        }
        //Ordering is very inefficient here. Should do all species for a single p2. Requires reordering.
//        for(s_ndx=0; s_ndx<system->num_stoch_species; s_ndx++){
//            off_diag_sum = 0.0; // keep track of the total off diagonal sum
//            for(n2=p1->neighbors->head; n2!=NULL; n2=n2->next){
//                p2 = n2->data;
//                diff_const = system->subdomain_diffusion_matrix[s_ndx*system->num_subdomains + (p2->type-1)];
//                //
//                //dist2 = particle_dist_sqrd(p1,p2);
//                // Eq (13-14), Drawert et al 2019
//                //ih = 1.0 / h;
//                //ihsq = ih * ih;
//                //wfd = h - sqrt(dist2);
//                //if(wfd <= 0.0){
//                //    continue; // outside support of basis function
//                //}
//                //wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih; //3D
//                // Eq 28 of Drawert et al 2019, Tartakovsky et. al., 2007, JCP
//                //D_i_j = -2.0*(p1->mass*p2->mass)/(p1->mass+p2->mass)*(p1->rho+p2->rho)/(p1->rho*p2->rho) * dist2 * wfd / (dist2+0.01*h*h);
//
//                if(diff_const > 0.0){
//                    //rdme->irD[irD_ndx++] = p2->id*system->num_stoch_species + s_ndx;
//                    //rdme->prD[prD_ndx++] = diff_const * D_i_j;  
//                    off_diag_sum += diff_const * n2->D_i_j;
//                }
//            }
//            //rdme->irD[irD_ndx++] = p1->id*system->num_stoch_species + s_ndx;
//            rdme->prD[prD_ndx++] = -1*off_diag_sum; 
//            //rdme->jcD[jcD_ndx++] = prD_ndx;
//        }
//    }
//    if(debug_flag){
//        printf("irD_ndx (%li) length rdme->irD (%li)\n",irD_ndx,irD_length*system->num_stoch_species);
//        printf("jcD_ndx (%li) length rdme->jcD (%li)\n",jcD_ndx,jcD_length*system->num_stoch_species);
//        printf("prD_ndx (%li) length rdme->prD (%li)\n",prD_ndx,prD_length*system->num_stoch_species);
//        if( prD_ndx != irD_ndx){
//            printf("Assembly: prD_ndx (%zu) != irD_ndx (%zu)\n",prD_ndx,irD_ndx);
//        }
//        if( irD_ndx != irD_length*system->num_stoch_species ){
//            printf("Assembly: irD_ndx (%zu) != irD_length*Mspecies (%li)\n", irD_ndx, irD_length*system->num_stoch_species);
//        }
//        char filename[256];
//        time_t seconds;
//        size_t i;
//        seconds = time(NULL);
//        sprintf(filename,"diffusion_matrix_%ld", seconds);
//        printf("Writing out diffusion matrix to '%s'\n",filename);
//        FILE*fp = fopen(filename,"w+");
//        fprintf(fp, "irD = [");
//        for(i=0;i<irD_ndx;i++){
//            if(i>0){ fprintf(fp,",");}
//            fprintf(fp, "%zu",rdme->irD[i]);
//        }
//        fprintf(fp, "]\n");
//        fprintf(fp, "jcD = [");
//        for(i=0;i<jcD_ndx;i++){
//            if(i>0){ fprintf(fp,",");}
//            fprintf(fp, "%zu",rdme->jcD[i]);
//        }
//        fprintf(fp, "]\n");
//        fprintf(fp, "prD = [");
//        for(i=0;i<prD_ndx;i++){
//            if(i>0){ fprintf(fp,",");}
//            fprintf(fp, "%e",rdme->prD[i]);
//        }
//        fprintf(fp, "]\n");
//        fprintf(fp, "D = scipy.sparse.csc_matrix(prD,irD,jcD)\n");
//        fclose(fp);
//
//    }

//    int num_subdomains;
//    const double*subdomain_diffusion_matrix;


//}

/**************************************************************************/
void nsm_core__destroy_diffusion_matrix(rdme_t*rdme){
    //free(rdme->irD);
    //free(rdme->jcD);
    //free(rdme->prD);
}


/**************************************************************************/
void nsm_core__take_step(system_t*system, double current_time, double step_size){
    rdme_t*rdme = system->rdme;
    double tt = current_time;
    double end_time = current_time + step_size;
    double totrate,cum,rdelta,rrdelta;
    int event,re,spec,errcode = 0;
    particle_t*subvol;
    size_t i,j = 0;
    size_t to_node,to_vol = 0;
    int dof,col;
    double old_rrate = 0.0,old_drate = 0.0;
    double rand1,rand2,cum2,old;
    double vol;

    /* Main loop. */
    while(tt <= end_time){

        /* Get the subvolume in which the next event occurred.
         This subvolume is on top of the heap. */
        //told = tt;
        //tt   = rdme->rtimes[0];
        //subvol = rdme->node[0];
        subvol = system->rdme->heap->head->data;
        tt = system->rdme->heap->head->tt;
        vol = (me->mass / me->rho);

        if(debug_flag){printf("nsm: tt=%e subvol=%i\n",tt,subvol->id);}
        /* First check if it is a reaction or a diffusion event. */
        totrate = subvol->rdme->srrate + subvol->rdme->sdrate;

//        if(totrate <= 0){ // Sanity check, is there a non-zero reaction and diffusion propensity
//            totrate = rdme->srrate[subvol]+rdme->sdrate[subvol];
//            if (totrate > 0.0)
//                rdme->rtimes[0] = -log(1.0-drand48())/totrate+tt;
//            else
//                rdme->rtimes[0] = INFINITY;
//            /* Update the heap. */
//            update(0,rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);
//            // go to the next element
//            continue;
//        }

        rand1 = drand48();

        if (rand1 <= subvol->rdme->srrate/totrate) { // use normalized floating point comparision
            /* Reaction event. */
            event = 0;

            /* a) Determine the reaction re that did occur (direct SSA). */
            double rand_rval = rand1 * subvol->rdme->srrate;
            for (re = 0, cum = subvol->rdme->rrate[0]; re < system->num_stoch_rxns && rand_rval > cum; re++, cum += subvol->rdme->rrate[re])
            ;
            if(re >= system->num_stoch_rxns){
                if(cum != subvol->rdme->srrate){
                    printf("Reaction propensity mismatch in voxel %i. re=%i, srrate[subvol]=%e cum=%e rand_rval=%e\n",subvol->id,re,subvol->rdme->srrate,cum,rand_rval);
                    rdelta = 0.0;
                    for (j=0;j<system->num_stoch_rxns; j++) {
                        rdelta += (subvol->rdme->rrate[j] = (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type));
                    }
                    subvol->rdme->srrate = rdelta;
                }
                if(subvol->rdme->srrate == 0.0){ continue; }


                double rand_rval2 = rand1 * subvol->rdme->srrate; // sum of propensitiess is not propensity sum, re-roll

                for (re = 0, cum = subvol->rdme->rrate[0]; re < system->num_stoch_rxns && rand_rval2 > cum; re++, cum += subvol->rdme->rrate[re])
                ;
                if(re >= system->num_stoch_rxns){ // failed twice, problems!
                    printf("Propensity sum overflow, rand=%e rand_rval=%e rand_rval2=%e srrate[%i]=%e cum=%e\n",rand1,rand_rval,rand_rval2,subvol->id,subvol->rdme->srrate,cum);
                    exit(1);
                }
            }
            if(debug_flag){printf("nsm: tt=%e subvol=%i sd=%i",tt,subvol->id,subvol->type);}
            if(debug_flag){printf("Rxn %i \n",re);}
            /* b) Update the state of the subvolume subvol and sdrate[subvol]. */
            for (i = rdme->jcN[re]; i < rdme->jcN[re+1]; i++) {
                int prev_val = subvol->xx[rdme->irN[i]];
                subvol->xx[rdme->irN[i]] += rdme->prN[i];
                if (subvol->xx[rdme->irN[i]] < 0){
                    errcode = 1;
                    printf("Netative state detected after reaction %i, subvol %i, species %zu at time %e (was %i now %i)\n",re,subvol->id,rdme->irN[i],tt,prev_val,subvol->xx[rdme->irN[i]]);

                    print_current_state(subvol,subvol->xx,system->num_stoch_species);
                    exit(1);
                }
                subvol->rdme->sdrate += subvol->rdme->Ddiag[rdme->irN[i]]*rdme->prN[i];
            }

            /* c) Recalculate srrate[subvol] using dependency graph. */
            for (i = rdme->jcG[system->num_stoch_species+re], rdelta = 0.0; i < rdme->jcG[system->num_stoch_species+re+1]; i++) {
                old = subvol->rdme->rrate[rdme->irG[i]];
                j = rdme->irG[i];
                rdelta +=
                (subvol->rdme->rrate[j] =
                 //(*rdme->rfun[j])(&rdme->xx[subvol*system->num_stoch_species],tt,rdme->vol[subvol],&rdme->data[subvol*rdme->dsize],rdme->sd[subvol])

                 (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type)
                 )-old;

            }
            subvol->rdme->srrate += rdelta;

            rdme->total_reactions++; /* counter */
        }
        else {
            /* Diffusion event. */
            event = 1;

            /* a) Determine which species... */
            double diff_rand = rand1 * subvol->rdme->sdrate;

            //for (spec = 0, dof = subvol*system->num_stoch_species, cum = rdme->Ddiag[dof]*rdme->xx[dof];
            //     spec < system->num_stoch_species && diff_rand > cum;
            //     spec++, cum += rdme->Ddiag[dof+spec]*rdme->xx[dof+spec]);
            for (spec = 0, cum = subvol->rdme->Ddiag[spec]*rdme->xx[spec]; 
                 spec < system->num_stoch_species && diff_rand > cum;
                 spec++, cum += subvol->rdme->Ddiag[spec]*rdme->xx[spec]);
            if(spec >= system->num_stoch_species){
                //printf("Diffusion species overflow\n");
                // try again, 'cum' is a better estimate of the propensity sum
                if(cum != subvol->rdme->srrate){
                    printf("Diffusion propensity mismatch in voxel %i. spec=%i, sdrate[subvol]=%e cum=%e diff_rand=%e\n",subvol->id,spec,subvol->rdme->sdrate,cum,diff_rand);
                    rdelta = 0.0;
                    for(j = 0; j < system->num_stoch_species; j++){
                        rdelta += subvol->rdme->Ddiag[j]*subvol->xx[j];
                    }
                    subvol->rdme->sdrate = rdelta;
                }
                if(subvol->rdme->sdrate == 0.0){ continue; }

                diff_rand = cum *rand1;
                for (spec = 0, cum = subvol->rdme->Ddiag[spec]*subvol->xx[spec];
                     spec < system->num_stoch_species && diff_rand > cum;
                     spec++, cum += subvol->rdme->Ddiag[spec]*subvol->xx[spec]);
                if(spec >= system->num_stoch_species){
                    spec--;
                    while(subvol->xx[spec] <= 0){
                        spec--;
                        if(spec <=0){
                            printf("Error: diffusion event in voxel %i was selected, but no molecues to move\n",subvol);
                            print_current_state(subvol,subvol->xx,system->num_stoch_species);
                            exit(1);
                        }
                    }
                }
            }


            /* b) and then the direction of diffusion. */
            double r2 = drand48();
            rand2 = r2 * subvol->rdme->Ddiag[spec];

            /* Search for diffusion direction. */
//            for (i = rdme->jcD[col], cum2 = 0.0; i < rdme->jcD[col+1]; i++)
//                if (rdme->irD[i] != col && (cum2 += rdme->prD[i]) > rand2)
//                    break;
//
//            /* paranoia fix: */
//            // This paranoia fix creates errors if the final rate has a zero propensity.  It can cause negative populations.
//            if (i >= rdme->jcD[col+1]){
//                //printf("Diffusion direction overflow\n");
//                // try again, 'cum2' is a better estimate of propensity sum
//                rand2 = r2*cum2;
//                for (i = rdme->jcD[col], cum2 = 0.0; i < rdme->jcD[col+1]; i++)
//                    if (rdme->irD[i] != col && (cum2 += rdme->prD[i]) > rand2)
//                        break;
//                if (i >= rdme->jcD[col+1]){
//                    i--;
//                }
//            }
//
//            to_node = rdme->irD[i];
//            to_vol = to_node/system->num_stoch_species;

            particle_t*dest_subvol = NULL
            particle_t*p2;
            neighbor_node_t*nn;
            cum2 = 0.0;
            for(nn=p->neighbors->head; nn!=NULL; nn=nn->next){
                p2 = nn->data;
                diff_const = system->subdomain_diffusion_matrix[spec*system->num_subdomains + (p2->type-1)];
                cum2 += nn->D_i_j * diff_const
                if(cum2 > rand2){
                    dest_subvol = p2;
                    break;
                }
            }
            if(nn==NULL){ //Diffusion direction overflow
                // try again, 'cum2' is a better estimate of propensity sum
                rand2 = r2*cum2;
                cum2 = 0.0;
                for(nn=p->neighbors->head; nn!=NULL; nn=nn->next){
                    p2 = nn->data;
                    diff_const = system->subdomain_diffusion_matrix[spec*system->num_subdomains + (p2->type-1)];
                    cum2 += nn->D_i_j * diff_const
                    if(cum2 > rand2){
                        dest_subvol = p2;
                        break;
                    }
                }
                if(nn==NULL){
                    printf("Error: overflow in trying to determine which destination voxel subvol=%i\n",subvol->id);
                    print_current_state(subvol,subvol->xx,system->num_stoch_species);
                    exit(1);
                }
            }


            /* c) Execute the diffusion event (check for negative elements). */
            subvol->xx[spec]--;
            if (subvol->xx[spec] < 0){
                    errcode = 1;
                    printf("Negative state detected after diffusion, voxel %i -> %zu, species %i at time %e\n",subvol,to_node,spec,tt);
                    print_current_state(subvol,rdme->xx,system->num_stoch_species);
                    exit(1);
            }

            dest_subvol->xx[spec]++;


            if(debug_flag){printf("nsm: tt=%e subvol=%i sd=%i",tt,subvol->id,subvol->type);}
            if(debug_flag){printf("Diff %i->%li\n",subvol->id,dest_subvol->id);}

            /* Save reaction and diffusion rates. */
            old_rrate = dest_subvol->rdme->srrate;
            old_drate = dest_subvol->rdme->sdrate;
            /* Recalculate the reaction rates using dependency graph G. */
            if (system->num_stoch_rxn > 0){
                for (i = rdme->jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < rdme->jcG[spec+1]; i++) {

// herehere
                    j = rdme->irG[i];
                    // update this voxel's reactions
                    old = subvol->rdme->rrate[j];
                    rdelta += 
                      (subvol->rdme->rrate[j] =
                        (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type)
                      )-old;


                    // update dest voxel's reactions
                    old = dest_subvol->rdme->rrate[j];
                    rrdelta += (dest_subvol->rdme->rrate[j] =
                        (*system->stoch_rxn_propensity_functions[j])(dest_subvol->xx,tt,vol,dest_subvol->data_fn,dest_subvol->type)
                      )-old;
                }

                subvol->rdme->srrate += rdelta;
                dest_subvol->rdme->srrate += rrdelta;
            }

            /* Adjust diffusion rates. */
            subvol->rdme->sdrate -= subvol->rdme->Ddiag[spec];
            dest_subvol->rdme->sdrate += dest_subvol->rdme->Ddiag[spec];

            rdme->total_diffusion++; /* counter */

        }

        /* Compute time to new event for this subvolume. */
        totrate = subvol->rdme->srrate+subvol->rdme->sdrate;
        if(totrate > 0.0){
            system->rdme->heap->head->tt = -log(1.0-drand48())/totrate+tt;
        }else{
            system->rdme->heap->head->tt = INFINITY;
        }
        /* Update the heap. */
        //update(0,rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);

        /* If it was a diffusion event, also update the other affected
         node. */
        if(event) {
            totrate = dest_subvol->rdme->srrate+dest_subvol->rdme->sdrate;
            if(totrate > 0.0) {
                if(!isinf(nn->tt)){
                    nn->tt =
                    (old_rrate+old_drate)/totrate*(nn->tt - tt)+tt;
                }else{
                    /* generate a new waiting time */
                    nn->tt = -log(1.0-drand48())/totrate+tt;
                }
            }else{
                nn->tt = INFINITY;
            }

            //update(rdme->heap[to_vol],rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);
        }

        // re-sort the heap
        ordered_list_sort(system->rdme->heap);
        //herehere

        /* Check for error codes. */
        if (errcode) {
            /* Cannot continue. Clear this solution and exit. */
            printf("Exiting due to errcode %i\n",errcode);
            print_current_state(subvol, rdme->xx,system->num_stoch_species);
            exit(1);
        }
    }



}

/**************************************************************************/
