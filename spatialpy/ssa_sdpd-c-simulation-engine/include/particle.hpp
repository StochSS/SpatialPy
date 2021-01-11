/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef particle_hpp
#define particle_hpp

#include <cstdlib>
#include <vector>
#include <queue>
#include "propensities.hpp"
// Include ANN KD Tree
#include <ANN/ANN.h>

namespace Spatialpy{

    struct Particle ;
    struct ParticleSystem ;
    struct NeighborNode ;
    struct EventNode ;

    struct Particle{
        Particle(unsigned int id=0) ;
	    unsigned int id;
	    int type;
	    double x[3];
	    double v[3];
	    double vt[3];
	    double mass;
	    double rho;
	    double nu;
	    int solidTag;
	    double bvf_phi;
	    double normal[3];
	    double F[3];
	    double Frho;
	    double Fbp[3];
	    // Data Function
	    double* data_fn;
	    // chem_rxn_system
	    unsigned int*xx; // populaion of discrete/stochastic species
	    double *C;  // concentration of chem species
	    double *Q;  // flux of chem species
	    // below here for simulation
	    std::vector<NeighborNode> neighbors ;

        // Moved from rdme_voxel_t
        double srrate;
        double* rrate;
        double sdrate;
        double* Ddiag;
        EventNode*heap_index;

	    double particle_dist(Particle p2);
	    double particle_dist_sqrd(Particle p2);
	    int add_to_neighbor_list(Particle neighbor, ParticleSystem system, double r2) ;

        // KD TREE FUNCTIONS
        void get_k_cleanup(ANNidxArray nn_idx, ANNdistArray dists) ;
        void search_cleanup(ANNpoint queryPt, ANNidxArray nn_idx, ANNdistArray dists) ;
        int get_k__approx(ParticleSystem system) ;
        int get_k__exact(ANNpoint queryPt, ANNdist dist, ParticleSystem system) ;
        void find_neighbors(ParticleSystem *system, bool use_exact_k=true) ;

            bool operator<(const Particle& p2){ 
                return x[0] > p2.x[0] ; 
            } 
    };

    struct EventNode{
    	Particle* data ;
    	double tt ;
        bool operator< (const EventNode& e2){ 
            return tt > e2.tt ; 
        } 
    };

    struct NeighborNode{
        Particle *data ;
        double dist ;
        double dWdr ;
        double D_i_j ;
	    NeighborNode(Particle *data, double dist, double dWdr, double D_i_j) ;

        bool operator<(const NeighborNode& n2){ 
            return dist > n2.dist ; 
        } 
    };

/**
    struct CompareTT { 
        bool operator()(EventNode const& e1, EventNode const& e2) 
        { 
            // return "false" if "e1" is ordered  
            // before "e2", for example: 
            return e1.tt > e2.tt ; 
        } 
    };
 
    struct CompareDist { 
        bool operator()(NeighborNode const& n1, NeighborNode const& n2) 
        { 
            // return "false" if "n1" is ordered  
            // before "n2", for example: 
            return n1.dist > n2.dist ; 
        } 
    };
 
    struct CompareX { 
        bool operator()(Particle const& p1, Particle const& p2) 
        { 
            // return "false" if "p1" is ordered  
            // before "p2", for example: 
            return p1.x[0] > p2.x[0] ; 
        } 
    };
**/
 
    struct ParticleSystem{
        ParticleSystem(size_t num_types, size_t num_chem_species, size_t num_chem_rxns, 
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn) ;
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
	std::vector<Particle> particles;
	std::vector<EventNode> event_v;
    // std::priority_queue<Particle> x_index;
    //std::priority_queue<EventNode> event_q;

        // Moved from rdme_t
        const size_t *irN;
        const size_t *jcN;
        const int *prN;
        const size_t *irG;
        const size_t *jcG;
        int initialized;
        long int total_reactions;
        long int total_diffusion;

	char boundary_conditions[3];
	bool static_domain;
	size_t num_types;
	size_t num_data_fn;

	ChemRxnFun* chem_rxn_rhs_functions;
	size_t num_chem_species;
	size_t num_chem_rxns;

	PropensityFun* stoch_rxn_propensity_functions;
	size_t num_stoch_species;
	size_t num_stoch_rxns;

	const char * const* species_names; 

	const double *subdomain_diffusion_matrix;
	//int *stochic_matrix;
	int *stoichiometric_matrix;
	double* gravity;

	void add_particle(Particle me);

        ANNkd_tree kdTree;
        ANNpointArray kdTree_pts;
        bool kdTree_initialized;
    };

}

/**
system_t* create_system(size_t num_types, size_t num_chem_species, size_t num_chem_rxns, 
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn);
particle_t* create_particle(int id);
**/

//bond_t* create_bond(particle p1, particle p2, double k, double rest_distance);


// Global flags
extern int debug_flag;


#endif //particle_h

