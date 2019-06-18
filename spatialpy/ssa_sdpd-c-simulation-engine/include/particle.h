/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef particle_h
#define particle_h

typedef struct __particle_t particle_t;
typedef struct __system_t system_t;
#include "linked_list.h"
#include "simulate_rdme.h"

struct __particle_t {
    unsigned int id;
    double x[3];
    double v[3];
    int type;
    double mass;
    double rho;
    double nu;
    int solidTag;
    // below here for simulation
    node* x_index;
    linked_list*neighbors;
};

struct __system_t {
    int dimension;
    double dt;
    unsigned int nt; 
    unsigned int current_step;
    unsigned int output_freq;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double h;
    double c0;
    double rho0;
    double P0;
    linked_list* particle_list;
    linked_list* x_index;
    char boundary_conditions[3];
    rdme_t*rdme;
    int static_domain;
};


linked_list* find_neighbors(particle_t* me, system_t* system);
system_t* create_system();
particle_t* create_particle(int id);
void add_particle(particle_t* me, system_t* system);
double particle_dist(particle_t* p1, particle_t*p2);
double particle_dist_sqrd(particle_t* p1, particle_t*p2);


// Global flags
extern int debug_flag;


#endif //particle_h

