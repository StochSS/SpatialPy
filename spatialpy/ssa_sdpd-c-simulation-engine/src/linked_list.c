/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#include <string.h>  // for strcmp and strcpy
#include <stdlib.h>  // for malloc and free
#include <stdio.h>   // for printf
#include "linked_list.h"
#include "particle.h"
#include <math.h>
#include "dSFMT/dSFMT.h"


//#define DEBUG_UPDATE


//constructor
linked_list_t* create_linked_list(){
    linked_list_t* ll = (linked_list_t*) malloc( sizeof(linked_list_t));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}
neighbor_list_t* create_neighbor_list(){
    neighbor_list_t* ll = (neighbor_list_t*) malloc( sizeof(neighbor_list_t));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}
ordered_list_t* create_ordered_list(){
    ordered_list_t* ll = (ordered_list_t*) malloc( sizeof(ordered_list_t));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}

// destructor
void destroy_linked_list( linked_list_t* ll ){
    // empty the linked list
    empty_linked_list(ll);
    // un-allocate the memory
    free(ll);
}
void destroy_neighbor_list( neighbor_list_t* ll ){
    // empty the linked list
    empty_neighbor_list(ll);
    // un-allocate the memory
    free(ll);
}
void destroy_ordered_list( ordered_list_t* ll ){
    // empty the linked list
    empty_ordered_list(ll);
    // un-allocate the memory
    free(ll);
}
// empty
void empty_linked_list( linked_list_t*ll){
    while( ll->count > 0){
        linked_list_delete( ll, ll->head );
    }
}
void empty_neighbor_list( neighbor_list_t*ll){
    while( ll->count > 0){
        neighbor_list_delete( ll, ll->head );
    }
}
void empty_ordered_list( ordered_list_t*ll){
    while( ll->count > 0){
        ordered_list_delete( ll, ll->head );
    }
}

// add a new node to the end of the linked list
node_t* linked_list_add( linked_list_t* ll, particle_t* data_in){
    node_t* n = (node_t *) malloc( sizeof(node_t) );
    n->data = data_in;
    n->next = NULL;
    n->prev = NULL;
    // Traverse the list to find the end node
    if(ll->head == NULL){
        ll->head = n;
        ll->tail = n;
    }else{
        ll->tail->next = n;
        n->prev = ll->tail;
        ll->tail = n;
    }
    // increase the size of the list
    ll->count++;
    return n;
}
neighbor_node_t* neighbor_list_add( neighbor_list_t* ll, particle_t* data_in){
    neighbor_node_t* n = (neighbor_node_t *) malloc( sizeof(neighbor_node_t) );
    n->data = data_in;
    n->next = NULL;
    n->prev = NULL;
    // Traverse the list to find the end node
    if(ll->head == NULL){
        ll->head = n;
        ll->tail = n;
    }else{
        ll->tail->next = n;
        n->prev = ll->tail;
        ll->tail = n;
    }
    // increase the size of the list
    ll->count++;
    return n;
}
ordered_node_t* ordered_list_add( ordered_list_t* ll, particle_t* data_in){
    ordered_node_t* n = (ordered_node_t *) malloc( sizeof(ordered_node_t) );
    n->data = data_in;
    n->next = NULL;
    n->prev = NULL;
    // Traverse the list to find the end node
    if(ll->head == NULL){
        ll->head = n;
        ll->tail = n;
    }else{
        ll->tail->next = n;
        n->prev = ll->tail;
        ll->tail = n;
    }
    // increase the size of the list
    ll->count++;
    return n;
}

// Delete a node from the linked list
void linked_list_delete( linked_list_t* ll, node_t* to_delete){
    node_t* prev_node;
    if( ll->head == NULL){
        printf("Error, linked_list_delete() empty list\n");
        return;
    }else if( to_delete == ll->head ){
        ll->head = ll->head->next;
    }else{
        for( prev_node=ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
            if(prev_node->next == to_delete){
                break;
            }
        }
        if( prev_node->next == NULL){
            printf("Error, linked_list_delete(), could not find item in list\n");
            return;
        }
        prev_node->next = to_delete->next;  // connect the list
        if(prev_node->next != NULL){ prev_node->next->prev = prev_node; }
    }

    //free and reduce size
    ll->count--;
    free(to_delete);
}
void neighbor_list_delete( neighbor_list_t* ll, neighbor_node_t* to_delete){
    neighbor_node_t* prev_node;
    if( ll->head == NULL){
        printf("Error, neighbor_list_t_delete() empty list\n");
        return;
    }else if( to_delete == ll->head ){
        ll->head = ll->head->next;
    }else{
        for( prev_node=ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
            if(prev_node->next == to_delete){
                break;
            }
        }
        if( prev_node->next == NULL){
            printf("Error, neighbor_list_t_delete(), could not find item in list\n");
            return;
        }
        prev_node->next = to_delete->next;  // connect the list
        if(prev_node->next != NULL){ prev_node->next->prev = prev_node; }
    }

    //free and reduce size
    ll->count--;
    free(to_delete);
}
void ordered_list_delete( ordered_list_t* ll, ordered_node_t* to_delete){
    ordered_node_t* prev_node;
    if( ll->head == NULL){
        printf("Error, ordered_list_t_delete() empty list\n");
        return;
    }else if( to_delete == ll->head ){
        ll->head = ll->head->next;
    }else{
        for( prev_node=ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
            if(prev_node->next == to_delete){
                break;
            }
        }
        if( prev_node->next == NULL){
            printf("Error, ordered_list_t_delete(), could not find item in list\n");
            return;
        }
        prev_node->next = to_delete->next;  // connect the list
        if(prev_node->next != NULL){ prev_node->next->prev = prev_node; }
    }

    //free and reduce size
    ll->count--;
    free(to_delete);
}

static inline void linked_list_sort__swap(node_t* a, node_t* b){
    particle_t*tmp = a->data;
    a->data = b->data;
    b->data = tmp;
    a->data->x_index = a ;
    b->data->x_index = b ;
}

static inline void neighbor_list_sort__swap(neighbor_node_t* a, neighbor_node_t* b){
    particle_t*tmp = a->data;
    a->data = b->data;
    b->data = tmp;
    double t2;
    t2 = a->dist;
    a->dist = b->dist;
    b->dist = t2;
    t2 = a->dWdr;
    a->dWdr = b->dWdr;
    b->dWdr = t2;
    t2 = a->D_i_j;
    a->D_i_j = b->D_i_j;
    b->D_i_j = t2;
}

static inline node_t *lastNode(node_t *root){
    while (root && root->next)
        root = root->next;
    return root;
}


static inline void linked_list_sort__swap(node_t* a, node_t* b){
    particle_t*tmp = b->data;
    b->data = a->data;
    a->data = tmp;
    a->data->x_index = a;
    b->data->x_index = b;
}
static inline node_t* linked_list_sort__partition_avg(node_t* min, node_t* max){
    double partitionValue = (min->data->x[0] + max->data->x[0])/2.0;
    node_t*left = min;
    node_t*right = max;
    while(left!=right && right->next != left){
        // from left, search for an element that is > partitionValue
        while(left!=right && right->next != left && left->data->x[0] <= partitionValue){
            left = left->next;
        }
        // from right, search for an element that is < partitionValue
        while(left!=right && right->next != left && right->data->x[0] > partitionValue){
            right = right->prev;
        }
        // swap elements
        if(left!=right && right->next != left){
            linked_list_sort__swap(left,right);
        }
    }
    return right;
}
void linked_list_sort__quicksort(node_t* min, node_t* max){
    if(min==NULL||max==NULL||min==max||max->next==min) return;
    node_t*pivot = linked_list_sort__partition_avg(min,max);
    linked_list_sort__quicksort(min, pivot->prev );
    linked_list_sort__quicksort(pivot, max );
}

void linked_list_sort(linked_list_t*ll){
    linked_list_sort__quicksort(ll->head, ll->tail);
    node_t*node = ll->head->next;
    while(node != NULL) {
        if(node->prev->data->x[0] > node->data->x[0]) {
            printf("===================================\n");
            printf("The list was not properly sorted\n");
            printf("===================================\n");
        }
        node = node->next;
    }
}


static inline void neighbor_list_sort__swap(neighbor_node_t* a, neighbor_node_t* b){
    particle_t*tmp = b->data;
    b->data = a->data;
    a->data = tmp;
    double t2;
    t2= b->dist;
    b->dist = a->dist;
    a->dist = t2;
    t2= b->dWdr;
    b->dWdr = a->dWdr;
    a->dWdr = t2;
    t2= b->D_i_j;
    b->D_i_j = a->D_i_j;
    a->D_i_j = t2;
}
static inline neighbor_node_t* neighbor_list_sort__partition_avg(neighbor_node_t* min, neighbor_node_t* max){
    double partitionValue = (min->dist + max->dist)/2.0;
    neighbor_node_t*left = min;
    neighbor_node_t*right = max;
    while(left!=right && right->next != left){
        // from left, search for an element that is > partitionValue
        while(left!=right && right->next != left && left->dist <= partitionValue){
            left = left->next;
        }
        // from right, search for an element that is < partitionValue
        while(left!=right && right->next != left && right->dist > partitionValue){
            right = right->prev;
        }
        // swap elements
        if(left!=right && right->next != left){
            neighbor_list_sort__swap(left,right);
        }
    }
    return right;
}
void neighbor_list_sort__quicksort(neighbor_node_t* min, neighbor_node_t* max){
    if(min==NULL||max==NULL||min==max||max->next==min) return;
    neighbor_node_t*pivot = neighbor_list_sort__partition_avg(min,max);
    neighbor_list_sort__quicksort(min, pivot->prev );
    neighbor_list_sort__quicksort(pivot, max );
}

void neighbor_list_sort(neighbor_list_t*ll){
    neighbor_list_sort__quicksort(ll->head, ll->tail->next);
    neighbor_node_t*node = ll->head->next;
    while(node != NULL) {
        if(node->prev->dist > node->dist) {
            printf("===================================\n");
            printf("The list was not properly sorted\n");
            printf("===================================\n");
        }
        node = node->next;
    }
}

static inline void ordered_list_sort__swap(ordered_node_t* a, ordered_node_t* b){
    particle_t*tmp = b->data;
    b->data = a->data;
    a->data = tmp;
    a->data->heap_index = a;
    b->data->heap_index = b;
    double t2;
    t2= b->tt;
    b->tt = a->tt;
    a->tt = t2;
}
static inline ordered_node_t* ordered_list_sort__partition_avg(ordered_node_t* min, ordered_node_t* max){
    double partitionValue = (min->tt + max->tt)/2.0;
    ordered_node_t*left = min;
    ordered_node_t*right = max;
    while(left!=right && right->next != left){
        // from left, search for an element that is > partitionValue
        while(left!=right && right->next != left && left->tt <= partitionValue){
            left = left->next;
        }
        // from right, search for an element that is < partitionValue
        while(left!=right && right->next != left && right->tt > partitionValue){
            right = right->prev;
        }
        // swap elements
        if(left!=right && right->next != left){
            ordered_list_sort__swap(left,right);
        }
    }
    return right;
}
void ordered_list_sort__quicksort(ordered_node_t* min, ordered_node_t* max){
    if(min==NULL||max==NULL||min==max||max->next==min) return;
    ordered_node_t*pivot = ordered_list_sort__partition_avg(min,max);
    ordered_list_sort__quicksort(min, pivot->prev );
    ordered_list_sort__quicksort(pivot, max );
}

void ordered_list_sort(ordered_list_t*ll){
    ordered_list_sort__quicksort(ll->head, ll->tail);
    ordered_node_t*node = ll->head->next;
    // while(node != NULL) {
    //     if(node->prev->tt > node->tt) {
    //         printf("===================================\n");
    //         printf("The list was not properly sorted\n");
    //         printf("===================================\n");
    //     }
    //     node = node->next;
    // }
}


// move a single element in an otherwise sorted list
void ordered_list_bubble_up_down(ordered_list_t*ll, ordered_node_t*n){
    ordered_node_t*n1 = n->next;

    // Remove node from current position
   
    ordered_node_t*before = n->prev;
    ordered_node_t*after  = n->next;
    if(before != NULL){ // If node before, connect that to one after
        before->next = after;
    }else{
        ll->head = ll->head->next; //else node is head
        ll->head->prev = NULL ;
    }
    if(after != NULL){ // If node after, connect that to one before
        after->prev = before; // else node is tail
    }else{ // if nothing after (is tail), set tail to node before
        ll->tail = ll->tail->prev;
        ll->tail->next = NULL ;
    }
    // Find new position
    //      if tt==inf, move to end
    if(isinf(n->tt) || n->tt >= ll->tail->tt){
        // move to end
        if(n == ll->head){
            ll->head = n->next ;
        }
        ll->tail->next = n;
        n->prev = ll->tail;
        n->next = NULL;
        ll->tail = n;
    }
    //  check to move to beginning
    else if(n->tt <= ll->head->tt || isinf(ll->head->tt)){
        // move to beginning
        ll->head->prev = n;
        n->next = ll->head;
        n->prev = NULL;
        ll->head = n;
    }
    //      Check if we move down, move down linearly
    else if(n->next != NULL && n->next->tt < n->tt){
        while(n1!=NULL && n1->tt < n->tt){
            n1=n1->next ;
        }
        n->next = n1;
        n->prev = n1->prev;
        n1->prev->next = n;
        n1->prev = n;
    }
    else {
        n1 = n->prev ;
        while(n1!=NULL && n1->tt > n->tt){
            n1=n1->prev ;
        }
        n->prev = n1;
        n->next = n1->next;
        n1->next->prev = n;
        n1->next = n;
    }


}

