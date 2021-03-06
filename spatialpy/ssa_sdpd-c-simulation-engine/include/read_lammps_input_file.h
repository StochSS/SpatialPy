/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef read_lammps_input_file_h
#define read_lammps_input_file_h
#include "linked_list.h"
#include "particle.h"

void read_lammps_input_file(const char*filename, system_t*system);


#endif //read_lammps_input_file_h

