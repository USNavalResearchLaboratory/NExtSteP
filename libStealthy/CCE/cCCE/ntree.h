#pragma once
#include "linkedList.h"
#include "arenaAlloc.h"

/*
    This structure defines a node in the CCE tree. Note that the
    children arrays are the first elements of the structure, this
    is done so that the offset to access them is the same reguardless
    of what other items are in the structs.  We abuse when
    implmenting sequence insertion and CCE calculation by casting the
    tRoot* to a tNode* because we only need to access the children.
*/
typedef struct treeNode {
    struct treeNode** children; // this must be the first element.
    unsigned int count;
} tNode;

typedef struct treeRoot {
    tNode** children; // this must be the first element.
    unsigned int branchingFactor;
    llRoot* layerWidth;
    arenaManager* arenaManagement;
} tRoot;

/*
    Creates a node for the tree and set all values to 0 or NULL;
    Input Args:
        int branchingFactor: the branching factor of the CCE tree
    Returns:
        tNode*: a pointer to the newly alocated node.
*/
tNode* createNode(tRoot* r, unsigned int branchingFactor);

/*
    Create a CCE tree with a given branching factor.
    Input Args:
        int branchinFactor: the branching factor of the tree
    retruns:
        tNode*: a pointer to the root element of the CCE Tree
*/
tRoot* createTree(unsigned int branchingFactor);

/*
    A wrapper around freeSubtree for freeing the entire tree
    Input Args:
        tNode* root: the root of the tree
*/
void freeTree(tRoot* root);