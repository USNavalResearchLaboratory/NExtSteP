/*
    An implmentation of a linked list.

*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "linkedList.h"


/*
    Creates the root of a linked list.  The root contains a count of
    length of the list and pointers to the front and back of the
    list.

    Returns:
        llRoot*: the ponter to the created root.
*/
llRoot* createLinkedList(void) {
    llRoot* r = (llRoot*)calloc(1, sizeof(llRoot));
    
    assert(r != NULL);
    
    return r;
}

/*
    Create a node of the linked list. A node contains a void* to some
    data and a pointer to the next item in the list.

    Input Args:
        void* value: The pointer you want inserted into the list.
    Returns:
        llNode*: The pointer to the created node.
*/
llNode* createllNode(void* value) {
    llNode* n = (llNode*)malloc(sizeof(llNode));
    
    assert(n != NULL);
    
    n->contents = value;
    n->next = NULL;
    return n;
}

/*
    Checks to see if the linked list is empty.

    Input Aargs:
        llRoot* r: A pointer to the root of the list to check.
    Returns:
        int: 1 if true, 0 if false.
*/
int isEmpty(llRoot* r) {
    assert(r != NULL);
    
    if(r->length) {
        return 0;
    } else {
        return 1;
    }
}

/*
    Insert node at the back of a list.

    Input Args:
        llRoot* r: The root of the linked list.
        void* value: The value to be inserted into the list.
*/
void insertllBack(llRoot* r, void* value) {
    llNode* n;
    
    assert(r != NULL);
    
    n = createllNode(value);
    if(isEmpty(r)) {
        r->head = n;
        r->tail = r->head;
    } else {
        r->tail->next = n;
        r->tail = n;
    }
    r->length += 1;
}

/*
    Insert node at the front of a list.

    Input Args:
        llRoot* r: The root of the linked list.
        void* value: The value to be inserted into the list.
*/
void insertllFront(llRoot* r, void* value) {
    llNode* n;
    
    assert(r != NULL);
    
    n = createllNode(value);
    if(isEmpty(r)) {
        r->head = n;
        r->tail = r->head;
    } else {
        n->next = r->head;
        r->head = n;
    }
    r->length += 1;
}

/*
    Frees the linked list. Only frees the linked list, anything
    pointed to by the list is not freed.

    Input Args:
        llRoot* r: the root of the linked list.
*/
void freeLL(llRoot* r) {
    llNode* n;
    llNode* tmp;
    
    assert(r != NULL);
    
    n = r->head;
    // walk the list and free each node
    while(n != NULL) {
        tmp = n;
        n = n->next;
        free(tmp);
    }
    free(r);
}

/*
    Walks the list and passes the pointer to llNode->contents to a
    user supplied void function which takes pointer to a void* 
    (void**) an argument.

    Input Args:
        llRoot* r: The root of the linked list.
        void (*func)(): The 
*/
void applyToList(llRoot* r, void (*func)(void**)) {
    llNode* n;    
    
    assert(r != NULL);
    assert(func != NULL);
    
    for(n = r->head; n != NULL; n = n->next) {
        // dereference and call the function pointer and pass the
        // pointer to contents to it.
        (*func)(&(n->contents));
    }
}

/*
    A generic print function for the linked list. Uses a user defined
    function to print each item of the list.  The user defined
    function takes a pointer to the contents of the node (void**) of
    the linked list. If the pointer to the user defined funciton is
    NULL the pointer to each element of the list is printed instead.

    A newline is printed at the end of the list.

    If the list is empty the word "empty" followed by a newline char
    will be printed.

    Input Args:
        llRoot* r: The root of the linked list.
        void (*printFunction)(): A pointer to a user defined print
            function.

*/
void printList(llRoot* r, void (*printFunction)(void**)) {
    llNode* n;

    assert(r != NULL);    
    
    if(isEmpty(r)) {
        printf("empty\n");
        return;
    }
    if(printFunction) {
        applyToList(r, printFunction);
        printf("\n");
    } else {
        for(n = r->head; n != NULL; n = n->next) {
            printf("%p ", n);
        }
        printf("\n");
    }
}