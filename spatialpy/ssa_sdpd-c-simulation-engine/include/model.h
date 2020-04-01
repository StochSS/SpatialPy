/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef model_h
#define model_h
#include "linked_list.h"
#include "particle.h"



void filterDensity(particle_t* me, linked_list* neighbors, system_t* system);

void pairwiseForce(particle_t* me, linked_list* neighbors, system_t* system);

void chemRxnFlux(particle_t* me, linked_list* neighbors, system_t* system);

void computeBoundaryVolumeFraction(particle_t* me, linked_list* neighbors, system_t* system);

void applyBoundaryCondition(particle_t* me, system_t* system);

void enforceVelocity(particle_t* me, system_t* system);

#endif //model_h
