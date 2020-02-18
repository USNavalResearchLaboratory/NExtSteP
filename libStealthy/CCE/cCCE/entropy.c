#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "linkedList.h"
#include "entropy.h"



/*
    Calculate the entropy of a linked list containing a layer of a
    CCE tree.
    Values are symbol frequencies within a layer
    Input Args:
        llRoot* r: the root of the linked list
        unsigned long numSymbols: the number of symbols in the
            alphabet in the case of CCE this is the number of
            times we have seen patterns go to that level.
    Returns:
        double: the entropy of that layer
*/
double entropy_ll(llRoot* r, unsigned long numSymbols) {
    double x;
    double entropy = 0;
    double dNumSymbols;
    llNode* n;

    assert(r != NULL);

    // Avoid repating expensive cast to floating point type
    dNumSymbols = (double)numSymbols;

    for(n = r->head; n != NULL; n = n->next) {
        x = ((double)(unsigned long)(n->contents)) / dNumSymbols;
        entropy -= x * log(x);
    }
    return entropy;
}

/*
    Calculate the entropy of an array containing a layer of a CCE
    tree.
    Values are symbol frequencies within a layer
    Input Args:
        unsigned long* r: the pointer to the array
        unsigned long numSymbols: the number of symbols in the
            alphabet in the case of CCE this is the number of
            times we have seen patterns go to that level.
        unsigned long length: the length of the pattern
    Returns:
        double: the entropy of that layer
*/
double entropy_arr(unsigned int* r, unsigned long numSymbols, unsigned long length) {
    double x, dNumSymbols;
    double entropy = 0;
    unsigned long n;

    assert(r != NULL);

    // Avoid repeating expensive cast to floating point type
    dNumSymbols = (double)numSymbols;

    for(n = 0; n < length; ++n) {
        x = (double)r[n] / dNumSymbols;
        entropy -= x * log(x);
    }

    return entropy;
}

/*
    Calculate what percentage of patterns are unique at this layer of
    the CCE tree.
    Input Args:
        llRoot* r: a pointer to a linked list containing the counts
            of a layer of a CCE tree.
    returns:
        double: the percentage of unique patterns in the layer.
*/
double fractionOfUniquePatterns_ll(llRoot* r) {
    double percentUnique = 0.0;
    unsigned long uniquePatterns = 0;
    unsigned long numberOfSequences = 0;
    llNode* n;

    assert(r != NULL);

    // Here we walk the list incrmenting the uniquePatterns counter
    // whenever we find a node that has only been visited once.
    for(n = r->head; n != NULL; n = n->next) {
        numberOfSequences += (unsigned long)n->contents;
        if((unsigned long)n->contents == 1) {
            ++uniquePatterns;
        }
    }

    // after counting how many patterns are unique we devide it by
    // the total number of patterns which in this case is represented
    // by the length of the linked list
    percentUnique = (double)uniquePatterns / (double)numberOfSequences;

    return percentUnique;
}

/*
    Given an array containing a layer of a CCE tree claculate what
    percentage of patterns are unique at this layer ofthe tree.
    Input Args:
        llRoot* r: 
            A pointer to an array containing the counts of a layer of
            a CCE tree.
        unsigned long length:
            The length of the layer.
    returns:
        double: the percentage of unique patterns in the layer.
*/
double fractionOfUniquePatterns_arr(unsigned int* r, unsigned long length) {
    double percentUnique = 0.0;
    unsigned int uniquePatterns = 0;
    unsigned long numberOfSequences = 0;
    unsigned long n;

    assert(r != NULL);

    // Here we walk the array incrmenting the uniquePatterns counter
    // whenever we find a node that has only been visited once.
    for(n = 0; n < length; ++n) {
        numberOfSequences += r[n];
        if(r[n] == 1) {
            ++uniquePatterns;
        }
    }

    // after counting how many patterns are unique we devide it by
    // the total number of patterns which in this case is the length.
    percentUnique = (double)uniquePatterns / (double)numberOfSequences;

    return percentUnique;
}

// int main(int argc, char const *argv[])
// {
//     // unsigned int values[10] = {2, 4, 3, 0, 0, 4, 2, 0, 2, 3};
//     unsigned int values[4] = {6,2,1,1};

//     double uniq = fractionOfUniquePatterns_arr(values, 4);
//     printf("%lf\n", uniq);
//     double ent = entropy_arr(values, 10, 4);
//     printf("%lf\n", ent);    
//     return 0;
// }
