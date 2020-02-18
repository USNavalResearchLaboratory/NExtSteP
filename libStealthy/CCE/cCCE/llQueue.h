#ifndef LLQUEUE
#define LLQUEUE 1

#include "linkedList.h"

/*
    Conveience wrappers around the linked list types to abstract them
    into queues.
*/
typedef llRoot queue;
typedef llNode qNode;

/*
    A conveience wrapper around createLinkedList to abstract it into
    a queue.
    Returns:
        queue*: A pointer to the created queue.
*/
queue* createQueue(void);

/*
    A conveience wrapper around insertllBack to abstract it into a
    queue and maintan compatablity with existing code.
    Input Args:
        queue* q: A pointer to the root of the queue
        void* contents: A pointer to the item you would like to
            insert into the queue.
    Returns:
        int: The length of the queue.
*/
int enqueue(queue* q, void* contents);


/*
    The only "real" function in this whole file. Removes the first
    item in the queue and returns is contents.
    Input Args:
        queue* q: A pointer to the root of the queue.
    Returns:
        void*: The cntents of the queue node.
*/
void* dequeue(queue* q);

/*
    A conveience wrapper around printList.
    Input Args:
        queue* q: A pointer to the root of the queue you wsh to print.
        void (*printFunction)(): A pointer to the function to be used
            to print each member of the queue.
*/
void printQ(queue* q, void (*printFunction)(void**));
#endif