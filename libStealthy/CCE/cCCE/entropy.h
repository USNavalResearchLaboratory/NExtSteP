#ifndef ENTROPY
#define ENTROPY 1
#include "linkedList.h"
/*
    Claculate the entropy of a linked list containing a layer of a
    CCE tree.
    Input Args:
        llRoot* r: the root of the linked list
        unsigned long numSymbols: the number of symbols in the 
            alphabet in the case of CCE this is the number of 
            times we have seen patterns go to that level.
    Returns:
        double: the entropy of tha layer
*/
double entropy_ll(llRoot* r, unsigned long numSymbols);

/*
    Claculate the entropy of an array containing a layer of a CCE
    tree.
    Input Args:
        unsigned long* r: the pointer to the array
        unsigned long numSymbols: the number of symbols in the 
            alphabet in the case of CCE this is the number of 
            times we have seen patterns go to that level.
        unsigned long length: the length of the pattern
    Returns:
        double: the entropy of that layer
*/
double entropy_arr(unsigned int* r, unsigned long numSymbols, unsigned long length);

/*
    Given an linked-list containing a layer of a CCE tree claculate 
    what percentage of patterns are unique at this layer ofthe tree.
    Input Args:
        llRoot* r: 
            A pointer to an array containing the counts of a layer of
            a CCE tree.
        unsigned long length:
            The length of the layer.
    returns:
        double: the percentage of unique patterns in the layer.
*/
double fractionOfUniquePatterns_ll(llRoot* r);

/*
    Given an array containing a layer of a CCE tree claculate what
    percentage of patterns are unique at this layer ofthe tree.
    Input Args:
        llRoot* r: 
            A pointer to a linked list containing the counts of a
            layer of a CCE tree.
        unsigned long length:
            The length of the layer.
    returns:
        double: The percentage of unique patterns in the layer.
*/
double fractionOfUniquePatterns_arr(unsigned int* r, unsigned long length);

#endif