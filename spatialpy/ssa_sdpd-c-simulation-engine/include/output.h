/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef output_h
#define output_h
#include "particle.h"



void output_csv(system_t*system, int current_step);

void output_vtk__sync_step(system_t*system, int current_step);
void output_vtk__async_step();

#endif // output_h

