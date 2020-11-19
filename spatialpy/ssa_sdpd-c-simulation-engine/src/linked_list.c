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

// search for a node by it's data field
/*
node_t* linked_list_search( linked_list_t* ll, char* search_string ){
    node_t* n;
    for( n=ll->head; n != NULL; n = n->next ){
        if( strcmp( n->data, search_string) == 0  ){
            break;
        }
    }
    if( n == NULL){
        return NULL;
    }
    // success, found the element
    return n;
}*/

// get node by index
/*
node_t* linked_list_get( linked_list_t* ll, int index){
    int count = 0;
    node_t* n = ll->head;
    if( ll->head == NULL){
        printf("Error, linked_list_get() empty list\n");
        return NULL;
    }
    while( count < index ){
        if(n->next == NULL){
            printf("Error, linked_list_get() list shorter than %i \n", index);
            return NULL;
        }
        n = n->next;
        count++;
    }
    return n;

}
*/

// remove and return first node on list
/*
node * linked_list_pop( linked_list * ll){
    node_t*n = ll->head;
    if( ll->head == NULL){
        return NULL;
    }
    ll->head = ll->head->next;
    ll->head->prev = NULL;
    ll->count--;
    return n;
}
*/


node_t* linked_list_sort__sub(node_t* head){
    node_t* min_node = head;
    node_t* before = NULL;
    node_t* ptr;
    node_t* tmp;
    if(head->next == NULL){
        return head;
    }
    for(ptr = head; ptr->next != NULL; ptr = ptr->next){
        if( ptr->next->data->x[0] < min_node->data->x[0] ){
            min_node = ptr->next;
            before = ptr;
        }
    }
    if( min_node != head ){
        tmp = head;
        head = min_node;
        before->next = min_node->next;
        if(min_node->next != NULL){ min_node->next->prev = before;}
        head->next = tmp;
        tmp->prev = head;
        head->prev = NULL;
    }
    head->next = linked_list_sort__sub(head->next);
    if(head->next != NULL){
        head->next->prev = head;
    }
    return head;
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
        while(left!=right && left->data->x[0] <= partitionValue){
            left = left->next;
        }
        // from right, search for an element that is < partitionValue
        while(left!=right && right->data->x[0] > partitionValue){
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

/*neighbor_node_t* neighbor_list_sort__sub(neighbor_node_t* head){
    neighbor_node_t* min_neighbor = head;
    neighbor_node_t* before = NULL;
    neighbor_node_t* ptr;
    neighbor_node_t* tmp;
    if(head->next == NULL){
        return head;
    }
    for(ptr = head; ptr->next != NULL; ptr = ptr->next){
        if( ptr->next->dist < min_neighbor->dist ){
            min_neighbor = ptr->next;
            before = ptr;
        }
    }
    if( min_neighbor != head ){
        tmp = head;
        head = min_neighbor;
        before->next = min_neighbor->next;
        if(min_neighbor->next != NULL){ min_neighbor->next->prev = before;}
        head->next = tmp;
        tmp->prev = head;
        head->prev = NULL;
    }
    head->next = neighbor_list_sort__sub(head->next);
    if(head->next != NULL){
        head->next->prev = head;
    }
    return head;
}*/

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
        while(left!=right && left->dist <= partitionValue){
            left = left->next;
        }
        // from right, search for an element that is < partitionValue
        while(left!=right && right->dist > partitionValue){
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
    // neighbor_node_t*node1 = ll->head;
    // while(node1 != NULL){
    //     printf("%f\n", node1->dist);
    //     node1 = node1->next;
    // }
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

/*ordered_node_t* ordered_list_sort__sub(ordered_node_t* head){
    ordered_node_t* min_ordered = head;
    ordered_node_t* before = NULL;
    ordered_node_t* ptr;
    ordered_node_t* tmp;
    if(head->next == NULL){
        return head;
    }
    for(ptr = head; ptr->next != NULL; ptr = ptr->next){
        if( ptr->next->tt < min_ordered->tt ){
            min_ordered = ptr->next;
            before = ptr;
        }
    }
    if( min_ordered != head ){
        tmp = head;
        head = min_ordered;
        before->next = min_ordered->next;
        if(min_ordered->next != NULL){ min_ordered->next->prev = before;}
        head->next = tmp;
        tmp->prev = head;
        head->prev = NULL;
    }
    head->next = ordered_list_sort__sub(head->next);
    if(head->next != NULL){
        head->next->prev = head;
    }
    return head;
}*/
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
        while(left!=right && left->tt <= partitionValue){
            left = left->next;
        }
        // from right, search for an element that is < partitionValue
        while(left!=right && right->tt > partitionValue){
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
    ordered_node_t*n1;
#ifdef DEBUG_UPDATE
    int cnt=0;
    for(n1=ll->head; n1!=NULL; n1=n1->next){ cnt++; }
    printf("before:");for(n1=ll->head; n1!=NULL; n1=n1->next){ printf("%i,",n1->data->id); }printf("\n");
    printf("ordered_list_bubble_up_down() id=%i tt=%e\tcnt=%i\n",n->data->id, n->tt,cnt);
    printf("\t");
#endif

    // Remove node from current position
    ordered_node_t*before = n->prev;
    ordered_node_t*after  = n->next;
    if(before != NULL){
        before->next = after;
#ifdef DEBUG_UPDATE
        printf("n->prev->tt=%e (id=%i)",n->prev->tt,n->prev->data->id);
#endif
    }else{
        ll->head = ll->head->next;
#ifdef DEBUG_UPDATE
        printf("n->prev=NULL ");
#endif
    }
    if(after != NULL){
        after->prev = before;
#ifdef DEBUG_UPDATE
        printf("n->next->tt=%e (id=%i)",n->next->tt,n->next->data->id);
#endif
    }else{
        ll->tail = ll->tail->prev;
#ifdef DEBUG_UPDATE
        printf("n->next=NULL ");
#endif
    }
#ifdef DEBUG_UPDATE
    printf("ll->head->tt=%e (id=%i)",ll->head->tt,ll->head->data->id);
    printf("ll->tail->tt=%e (id=%i)",ll->tail->tt,ll->tail->data->id);
    fflush(stdout);
    int mvcnt=0;

    printf("removed:");for(n1=ll->head; n1!=NULL; n1=n1->next){ printf("%i,",n1->data->id); }printf("\n");
#endif
    // Find new position
    //      if tt==inf, move to end
    if(isinf(n->tt) || n->tt >= ll->tail->tt){
        // move to end
        ll->tail->next = n;
        n->prev = ll->tail;
        n->next = NULL;
        ll->tail = n;
#ifdef DEBUG_UPDATE
        printf("\tmoved to end\n");
        int ecnt=0;
        for(n1=ll->head; n1!=NULL; n1=n1->next){ ecnt++; } 
        printf("after:");for(n1=ll->head; n1!=NULL; n1=n1->next){ printf("%i,",n1->data->id); }printf("\n");
        if(cnt!=ecnt){printf("count mismatch cnt=%i ecnt=%i\n",cnt,ecnt);exit(1);}
#endif
        return;
    }
    //  check to move to beginning
    else if(n->tt <= ll->head->tt || isinf(ll->head->tt)){
        // move to beginning
        ll->head->prev = n;
        n->next = ll->head;
        n->prev = NULL;
        ll->head = n;
#ifdef DEBUG_UPDATE
        printf("\tmoved to beginning\n");
        int ecnt=0;
        for(n1=ll->head; n1!=NULL; n1=n1->next){ ecnt++; } 
        if(cnt!=ecnt){printf("count mismatch cnt=%i ecnt=%i\n",cnt,ecnt);exit(1);}
#endif
        return;
    }
    //      Check if we move down, move down linearly
    else if(n->next != NULL && n->next->tt < n->tt){
        for(n1=n->next; n1!=NULL; n1=n1->next){ // find position before n1
#ifdef DEBUG_UPDATE
            mvcnt++;
#endif
            if(n1->next == NULL || n1->tt >= n->tt){
                n->next = n1;
                n->prev = n1->prev;
                n1->prev->next = n;
                n1->prev = n;
#ifdef DEBUG_UPDATE
                printf("\tmoved down %i\n", mvcnt);
                //printf("node: n->id=%i, n->tt=%e\n",n->data->id, n->tt);fflush(stdout);
                //printf("list:\n");fflush(stdout);
                //for(n1=ll->head;n1!=NULL;n1=n1->next){
                //    printf("\tid=%i tt=%e\n",n1->data->id, n1->tt);fflush(stdout);
                //}
                int ecnt=0;
                for(n1=ll->head; n1!=NULL; n1=n1->next){ ecnt++; } 
                if(cnt!=ecnt){printf("count mismatch cnt=%i ecnt=%i\n",cnt,ecnt);exit(1);}
#endif
                return;
            }
        }
    }
    //      check if we move up, move up linearly
    //else if(n->prev != NULL && n->prev->tt >= n->tt){
    else {
        for(n1=n->prev; n1!=NULL; n1=n1->prev){ // find position after n1 
            if(n1->prev == NULL || n1->tt <= n->tt){
                n->prev = n1;
                n->next = n1->next;
                n1->next->prev = n;
                n1->next = n;
#ifdef DEBUG_UPDATE
                printf("\tmoved up %i\n", mvcnt);
                int ecnt=0;
                for(n1=ll->head; n1!=NULL; n1=n1->next){ ecnt++; } 
                if(cnt!=ecnt){printf("count mismatch cnt=%i ecnt=%i\n",cnt,ecnt);exit(1);}
#endif
                return;
            }
#ifdef DEBUG_UPDATE
            mvcnt++;
#endif
        }
    }

    printf("ERROR, should not get here, ordered_list_bubble_up_down, node not inserted.\n");
#ifdef DEBUG_UPDATE
    printf("node: n->id=%i, n->tt=%e\n",n->data->id, n->tt);
    printf("list:\n");
    for(n1=ll->head;n1!=NULL;n1=n1->next){
        printf("\tid=%i tt=%e\n",n1->data->id, n1->tt);
    }
#endif
    exit(1);

}

