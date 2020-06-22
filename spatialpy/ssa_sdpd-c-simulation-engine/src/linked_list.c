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


//constructor
linked_list* create_linked_list(){
    linked_list* ll = (linked_list*) malloc( sizeof(linked_list));
    ll->count = 0;
    ll->head = NULL;
    ll->tail = NULL;
    return ll;
}

// destructor
void destroy_linked_list( linked_list* ll ){
    // empty the linked list
    empty_linked_list(ll);
    // un-allocate the memory
    free(ll);
}

void empty_linked_list( linked_list*ll){
    while( ll->count > 0){
        linked_list_delete( ll, ll->head );
    }
}

// add a new node to the end of the linked list
node* linked_list_add( linked_list* ll, particle_t* data_in){
    node* n = (node *) malloc( sizeof(node) );
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
void linked_list_delete( linked_list* ll, node* to_delete){
    node* prev_node;
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
    free_node(to_delete);
}

// search for a node by it's data field
/*
node* linked_list_search( linked_list* ll, char* search_string ){
    node* n;
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
node* linked_list_get( linked_list* ll, int index){
    int count = 0;
    node* n = ll->head;
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

// remove and return first node on list
node * linked_list_pop( linked_list * ll){
    node*n = ll->head;
    if( ll->head == NULL){
        return NULL;
    }
    ll->head = ll->head->next;
    ll->head->prev = NULL;
    ll->count--;
    return n;
}

void free_node(node*n){
    free(n);
}

node* linked_list_sort__sub(node* head, int sort_ndx){
    node* min_node = head;
    node* before = NULL;
    node* ptr;
    node* tmp;
    if(head->next == NULL){
        return head;
    }
    for(ptr = head; ptr->next != NULL; ptr = ptr->next){
        if( ptr->next->data->x[sort_ndx] < min_node->data->x[sort_ndx] ){
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
    head->next = linked_list_sort__sub(head->next,sort_ndx);
    if(head->next != NULL){
        head->next->prev = head;
    }
    return head;
}

void linked_list_sort(linked_list*ll, int sort_ndx){
    ll->head = linked_list_sort__sub(ll->head, sort_ndx);
}













