#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "ntree.h"
#include "linkedList.h"
#include "arenaAlloc.h"
#ifndef ARENA_SIZE
    #define ARENA_SIZE 2048
#endif

/*
    Creates a node for the tree and set all values to 0 or NULL;
    Input Args:
        int branchingFactor: the branching factor of the CCE tree
    Returns:
        tNode*: a pointer to the newly alocated node.
*/
tNode* createNode(tRoot* r, unsigned int branchingFactor) {
    unsigned long branchingSize;

    assert(branchingFactor > 0);

    branchingSize = branchingFactor * sizeof(tNode*);

    tNode* t = (tNode*)arenaAlloc(r->arenaManagement, sizeof(tNode));

    assert(t != NULL);

    t->count = 0;
    t->children = (tNode**)arenaAlloc(r->arenaManagement, branchingSize);
    // set all children to NULL
    memset(t->children, 0, branchingSize);
    return t;
}

/*
    Create a tree with a given branching factor.
    Input Args:
        int branchinFactor: the branching factor of the tree
    retruns:
        tNode*: a pointer to the root element of the CCE Tree
*/
tRoot* createTree(unsigned int branchingFactor) {
    tRoot* r = (tRoot*)malloc(sizeof(tRoot));

    assert(r != NULL);

    r->children = (tNode**)calloc(branchingFactor, sizeof(tNode*));
    r->branchingFactor = branchingFactor;
    r->layerWidth = createLinkedList();
    r->arenaManagement = setupArena(ARENA_SIZE);
    return r;
}

/*
    Frees the members of the tRoot struct.
    Input Args:
        tNode* root: the root of the tree
*/
void freeTree(tRoot* root) {

    assert(root != NULL);

    freeLL(root->layerWidth);
    teardownArena(root->arenaManagement);
    free(root->children);
    free(root);
}
