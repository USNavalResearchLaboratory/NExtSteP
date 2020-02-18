#ifndef INCLUDEDCCE
#define INCLUDEDCCE 1
#include "ntree.h"
#include "linkedList.h"

/*
    Insert a sequence into CCE tree.

    Input Args:
        tNode* root: the root of the CCE tree.
        int* sequence: a pointer to memory containg the binned sequence.
        int length: the length of the sequence.
*/
void insertSequence(tRoot* root, int* sequence, int length);

/* 
    Function used by the generic print list function to print each
    node in the layer.

    Input Args:
        void** n: The pointer to the contents of the linked list
            node. In this case it is a pointer to an unsigned long
            containg the count of how many times a node was visited.
*/
void printLayer(void** n);

/*
    Generate a list of Corrected Conditional Entropys for each layer
    of a given CCE tree.

    Input Args:
        tNode* root: the root of the CCE tree.
    Returns:
        double*: A pointer an array containing the calculated CCE
        values, this array is NAN terminated.
*/
double* calcCCEs(tRoot* root);

/*
    Fills memory pointed to by sequenceMemory with length random
    values between 0 and 4. Used for testing.

    Input Args:
        int length: The length of the memory to fill.
    Output Args:
        int* sequenceMemory: The ponter to the memory to be filled.
*/
void genRandomSequence(int* sequenceMemory, int length);
#endif