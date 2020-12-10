/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#include "particle.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <queue>

namespace Spatialpy{
    ParticleSystem::ParticleSystem(size_t num_types, size_t num_chem_species, size_t num_chem_rxns, 
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn){
        dimension = 3;
        boundary_conditions[0] = 'n';
        boundary_conditions[1] = 'n';
        boundary_conditions[2] = 'n';
        static_domain = 0;
        num_stoch_species = num_stoch_species;
        num_stoch_rxns = num_stoch_rxns;
        num_chem_species = num_chem_species;
        num_chem_rxns = num_chem_rxns;
        num_types = num_types;
        gravity = calloc(3,sizeof(double));
        num_data_fn = num_data_fn;
    }

    void ParticleSystem::add_particle(Particle me){
	x_index.push(me) ;
	particles.push_back(me) ;
	// TODO: does this need to change? 
	Q = (double*) calloc(system->num_chem_species, sizeof(double));
	C = (double*) calloc(system->num_chem_species, sizeof(double));
	data_fn = (double*) calloc(system->num_data_fn, sizeof(double));
    }

    Particle::Particle(int id){
	id = id;
	nu = 0.01; 
	mass = 1;
	rho = 1;
	solidTag = 0;
	x[0] = x[1] = x[2] = 0.0;
	v[0] = v[1] = v[2] = 0.0;
    }

    double Particle::particle_dist(Particle p2){
        double a = x[0] - p2.x[0];
        double b = x[1] - p2.x[1];
        double c = x[2] - p2.x[2];
        return sqrt( a*a + b*b + c*c);
    }

    double Particle::particle_dist_sqrd(Particle p2){
        double a = x[0] - p2.x[0];
        double b = x[1] - p2.x[1];
        double c = x[2] - p2.x[2];
        return ( a*a + b*b + c*c);
    }

    int Particle::add_to_neighbor_list(Particle neighbor, ParticleSystem system){
        double a = x[0] - neighbor.x[0];
	double b = x[1] - neighbor.x[1];
	double c = x[2] - neighbor.x[2];
	double r2 =  ( a*a + b*b + c*c);
	double r = sqrt(r2);

	if(r > system->h){ return 0; } // do not add, out side support radius

	// calculate dWdr
	double h = system.h;
	double R = r / h;
	double alpha = 105 / (16 * M_PI * h * h * h); // 3D
	double dWdr = alpha * (-12 * r / (h * h)) * ((1 - R)* (1 - R));
	// calculate D_i_j

	// Eq (13-14), Drawert et al 2019
	double ih = 1.0 / h;
	double ihsq = ih * ih;
	double dhr = h - r;
	double wfd = -25.066903536973515383e0 * dhr * dhr * ihsq * ihsq * ihsq * ih; //3D
	// Eq 28 of Drawert et al 2019, Tartakovsky et. al., 2007, JCP
	double D_i_j = -2.0*(mass*neighbor.mass)/(mass+neighbor.mass)*(rho+neighbor.rho)/(rho*neighbor.rho) * r2 * wfd / (r2+0.01*h*h);

	if(isnan(D_i_j)){
	    printf("Got NaN calculating D_i_j for me=%i, neighbor=%i\n",id, neighbor->id);
	    printf("r=%e ",r);
	    printf("h=%e ",h);
	    printf("alpha=%e ",alpha);
	    printf("dWdr=%e ",dWdr);
	    Particle p = me;
	    printf("mass=%e ",p->mass);
	    printf("rho=%e ",p->rho);
	    p = neighbor;
	    printf("n->mass=%e ",p->mass);
	    printf("n->rho=%e ",p->rho);

	    exit(1);
	}

	neighbors.push({neighbor, r, dWdr, D_i_j}) ;

	return 1;
    }

    void Particle::find_neighbors(ParticleSystem system){
	//clean out previous neighbors
	//printf("find_neighbors.empty_linked_list\n");
	//TODO: empty_neighbor_list(neighbors);
	// search for points forward: (assumes the list is sorted ascending)
	//printf("find_neighbors.search forward\n");
	for(Particle n : system.x_index){
	    if(n.data.x[0] > (x[0] + system.h)) break; //stop searching
	    if( (n.data.x[1] > (x[1] + system.h)) || (n.data.x[1] < (x[1] - system.h) ) ) continue;
		//TODO: Should these be breaks instead of continue?
	    if( (n.data.x[2] > (x[2] + system.h)) || (n.data.x[2] < (x[2] - system.h) ) ) continue;
	    add_to_neighbor_list(n, system);
	    if(debug_flag>2){ 
		printf("find_neighbors(%i) forward found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
		    id,n.data.id, particle_dist(n.data),
		    x[0] - n.data.x[0],
		    x[1] - n.data.x[1],
		    x[2] - n.data.x[2]);
	    }
	}

	//search for points backward
	for(std::vector<Particle>::iterator n = x_index.end(); n!=x_index.begin(); --n){
	    if(n->data->x[0] < (x[0] - system->h)){
		break; //stop searching backward
	    }
	    if( (n->data.x[1] > (x[1] + system.h)) || (n->data.x[1] < (x[1] - system.h) ) ){
		continue;
	    }
	    if( (n->data.x[2] > (x[2] + system.h)) || (n->data.x[2] < (x[2] - system.h) ) ){ 
		continue;
	    }
	    add_to_neighbor_list(n->data, system);
	    if(debug_flag>2){ 
		printf("find_neighbors(%i) backwards found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
		    id,n->data.id, particle_dist(n->data),
		    x[0] - n->data.x[0],
		    x[1] - n->data.x[1],
		    x[2] - n->data.x[2]);
	    }
	}
    }

}

