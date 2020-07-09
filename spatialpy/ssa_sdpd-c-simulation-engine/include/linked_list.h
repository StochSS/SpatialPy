/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef linked_list_h
#define linked_list_h
//include <pthread.h>
#include "pthread_barrier.h"

typedef struct linked_list_t linked_list;
typedef struct node_t node;
#include "particle.h"

// Create data structure for a node of the list
struct node_t {
    particle_t* data;
    node* next;
    node* prev;
};

struct neighbor_node_t {
    particle_t* data;
    double dist;
    double dWdr;
    double D_i_j;
    neighbor_node_t* next;
    neighbor_node_t* prev;
};

// Data structre for the linked list type

struct linked_list_t {
    node* head;
    node* tail;
    int count;
    //pthread_mutex_t mutex; // not used
};

// Functions to manipulate the linked list

//constructor
linked_list* create_linked_list();

void empty_linked_list( linked_list*ll);

// destructor
void destroy_linked_list( linked_list* ll );

// add a new node to the end of the linked list
node* linked_list_add( linked_list* ll, particle_t* data_in);

// Delete a node from the linked list
void linked_list_delete( linked_list* ll, node* to_delete);

// search for a node by it's data field
//node* linked_list_search( linked_list* ll, char* search_string );

// get node by index
node* linked_list_get( linked_list* ll, int index);

// get first node on list and remove it from list
node * linked_list_pop( linked_list * ll);

// free a node
void free_node(node*n);

// in-place bubble sort
void linked_list_sort(linked_list*ll, int sort_ndx);


#endif /* linked_list_h*/
