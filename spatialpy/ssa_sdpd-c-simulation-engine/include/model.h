/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef model_h
#define model_h
#include "linked_list.h"
#include "particle.h"



void filterDensity(particle_t* me, system_t* system);

void pairwiseForce(particle_t* me, system_t* system);

void computeBoundaryVolumeFraction(particle_t* me, system_t* system);

void applyBoundaryVolumeFraction(particle_t* me, system_t* system);


#endif //model_h
