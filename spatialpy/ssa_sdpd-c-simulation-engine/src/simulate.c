#include "linked_list.h"
#include "output.h"
#include "particle.h"
#include "simulate_rdme.h"
#include "model.h"
#include <errno.h>
#include "pthread_barrier.h"
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


// Step 1/3: First part of time step computation
void take_step1(particle_t* me, system_t* system, unsigned int step)
{
    int i;

    // Step 1.1: Enforce initial velocity conditions at step 0
    // This is now done directly via python interface
    //if (step == 0)
    //    enforceVelocity(me, system);

    // Step 1.2: Predictor step
    
    // Update half-state
    if (me->solidTag == 0 && system->static_domain == 0) {
       for (i = 0; i < 3; i++) {
           // Update velocity using forces
            me->v[i] = me->v[i] + 0.5 * system->dt * me->F[i];
            // Update transport velocity using background pressure force
            me->vt[i] = me->v[i] + 0.5 * system->dt * me->Fbp[i];
            // Update position using previous velocity
            me->x[i] = me->x[i] + system->dt * me->vt[i];
        }
        // Update density using continuity equation 
        me->rho = me->rho + 0.5 * system->dt * me->Frho;
    }
    // update half-state of chem rxn
    for(i=0; i< system->num_chem_species; i++){
        me->C[i] += me->Q[i] * system->dt * 0.5;
    }

    // Apply boundary conditions
    applyBoundaryConditions(me, system);



    // Step 1.3: Clean forces
    for (i = 0; i < 3; i++) {
        // Clean momentum force
        me -> F[i] = 0.0;
        // Clean background pressure force
        me -> Fbp[i] = 0.0;
    }
    // Clean mass flux term
    me -> Frho = 0.0;
    // Clean chem rxn flux
    for(i=0; i< system->num_chem_species; i++){
        me->Q[i] = 0.0;
    }

}

// Step 2/3: Update the 'F' field of each partile attached to a bond
// void compute_bond_forces(bond_t*bond, system_t*system, unsigned int step){
//}


// Step 2/3: Compute forces
void compute_forces(particle_t* me, system_t* system, unsigned int step)
{

    // Step 2.1: Build neighbor list at first step
    if (system->static_domain) {
        if (step == 0) {
            find_neighbors(me, system);
        }
        return;
    }

    // Step 2.2: Find nearest neighbors
    find_neighbors(me, system);

    // Step 2.3: Compute forces
    pairwiseForce(me, me->neighbors, system);
}


// Step 3/3: Compute the final state
void take_step2(particle_t* me, system_t* system, unsigned int step)
{
    int i;

    // Step 3.1: Corrector step
    if (me->solidTag == 0 && system->static_domain == 0) {
        for (i = 0; i < 3; i++) {
            // Update velocity using forces
            me->v[i] = me->v[i] + 0.5 * system->dt * me->F[i];
        }

        // Update density using continuity equation and Shepard filter 
        if (step % 20 == 0) {
            filterDensity(me, me->neighbors, system);
            me->rho = me->rho + 0.5 * system->dt * me->Frho;
        }
        else
            me->rho = me->rho + 0.5 * system->dt * me->Frho;
    }
    else {
        // Filter density field (for fixed solid particles)
        if (step % 20 == 0)
            filterDensity(me, me->neighbors, system);
    }


    // Step 3.2: Compute boundary volume fractions (bvf)
    // Step 3.3: Apply BVF 
    if (me->solidTag == 0) {
        computeBoundaryVolumeFraction(me, me->neighbors, system);
        applyBoundaryVolumeFraction(me, system);
    }

    //  Solve deterministic/stochastic reaction-diffusion system
    // update half-state of chem rxn
    for(i=0; i< system->num_chem_species; i++){
        me->C[i] += me->Q[i] * system->dt * 0.5;
    }
    // Apply boundary conditions
    applyBoundaryConditions(me, system);

}
