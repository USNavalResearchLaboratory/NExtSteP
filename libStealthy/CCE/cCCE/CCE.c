// #define DEBUG_CCE 1
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "linkedList.h"
#include "llQueue.h"
#include "entropy.h"
#include "ntree.h"
#include "CCE.h"


/*
    Insert a sequence into CCE tree.

    Input Args:
        tNode* root: the root of the CCE tree.
        int* sequence: a pointer to memory containg the binned sequence.
        int length: the length of the sequence.
*/
void insertSequence(tRoot* root, int* sequence, int length) {
    int i;
    int currentBin;
    tNode* currentNode;
    llNode* layerWidth;
    
    assert(root != NULL);
    assert(sequence != NULL);

    // Ensure we have enough nodes in the ll for the length of the input.
    if(root->layerWidth->length < length){
        for(i=root->layerWidth->length; i < length; ++i){
            insertllBack(root->layerWidth, (void*)0);
        }
    }

    currentNode = (tNode*)root;
    layerWidth = root->layerWidth->head;
    for(i=0; i<length; ++i) {
        currentBin = sequence[i];
        // if this path hasn't been taken before create a new node
        // and increase the layer width
        if(currentNode->children[currentBin] == NULL) {
            currentNode->children[currentBin] = createNode(root, root->branchingFactor);
            layerWidth->contents += 1;
        }
        // move to the child
        currentNode = currentNode->children[currentBin];
        currentNode->count += 1;
        // printf("bin %d, count %d\n", currentBin, currentNode->count);
        layerWidth = layerWidth->next;
    }
}

/*
    Function used by the generic print list function to print each
    node in the layer.

    Input Args:
        void** n: The pointer to the contents of the linked list
            node. In this case it is a pointer to an unsigned long
            containg the count of how many times a node was visited.
*/
void printLayer(void** n) {
    assert(n != NULL);
    printf("%lu ", *(unsigned long*)n);
}

/*
    Generate a list of Corrected Conditional Entropys for each layer
    of a given CCE tree.

    Input Args:
        tNode* root: the root of the CCE tree.
    Output Args:
        llRoot* cceLL: a linked list contaning the CCE vaues for each layer
*/
double* calcCCEs(tRoot* root) {
    int i, ii, depth;
    void* tmp;
    double* CCE;
    double conditionalEntropy;
    double previousLayerEntropy;
    double firstOrderEntropy;
    unsigned long rowSum;
    unsigned long layerWidth;
    unsigned long maxWidth;
    unsigned int* bfsLayer;
    double layerEntropy;
    double layerPercUnique;

    assert(root != NULL);

    tNode* currentNode = NULL;
    queue* currentLayer = createQueue();
    queue* nextLayer = createQueue();
    llNode* currentLayerWidth = root->layerWidth->head;
    layerEntropy = 0.0;
    layerPercUnique = 0.0;
    previousLayerEntropy = 0.0;
    depth = 0;
    maxWidth = (unsigned long)root->layerWidth->tail->contents;


    bfsLayer = (unsigned int*)malloc(sizeof(unsigned int) * maxWidth);
    CCE = (double*)malloc(sizeof(double) * root->layerWidth->length);

    assert(bfsLayer != NULL);
    assert(CCE != NULL);

    if(!enqueue(nextLayer, root)) {
        printf("unable to enqueue\n");
    }
    #ifdef DEBUG_CCE
        printf("CCE, Conditional Entropy, Correction Factor\n");
    #endif
    
    // here we perfrom a BFS of the tree and calculate the CCE for
    // each layer.
    for(depth=0; depth < root->layerWidth->length; ++depth) {
        layerWidth = (unsigned long)currentLayerWidth->contents;
        rowSum = 0;
        ii = 0;
        // swap queues
        tmp = (void*)currentLayer;
        currentLayer = nextLayer;
        nextLayer = (queue*)tmp;
        
        // Here we populate the current layer to be sent off to the
        // processing thread tree.  We also populate the nextLayer
        // with pointers to the next layer.
        while(!isEmpty(currentLayer)) {
            currentNode = (tNode*)dequeue(currentLayer);

            assert(currentNode != NULL);
            
            // check for child nodes, if they exist enqueue them.
            for(i = 0; i < root->branchingFactor; ++i) {
                tmp = (void*)currentNode->children[i];
                if(tmp != NULL) {
                    if(!enqueue(nextLayer, tmp)) {
                        printf("unable to enqueue\n");
                    }
                    bfsLayer[ii] = ((tNode*)tmp)->count;
                    rowSum += bfsLayer[ii];
                    ++ii;
                }
            }
        }

        // calculate CCE
        layerEntropy = entropy_arr(bfsLayer, rowSum, layerWidth);
        // set first order entropy
        if(depth == 0) {
            firstOrderEntropy = layerEntropy;
        }

        // calculate the conditional entropy
        conditionalEntropy = layerEntropy - previousLayerEntropy;
        previousLayerEntropy = layerEntropy;

        layerPercUnique = fractionOfUniquePatterns_arr(bfsLayer, layerWidth);

        CCE[depth] = conditionalEntropy + (firstOrderEntropy * layerPercUnique);

        #ifdef DEBUG_CCE
            printf("%f, %f, %f\n", CCE[depth], conditionalEntropy, (firstOrderEntropy * layerPercUnique));
        #endif

        // If all nodes on a layer are unique and we don't have any larger
        // layers below us we can be sure that of all layers below us have the
        // same CCE.
        if(layerPercUnique == 1.0){
            break;
        }
        // setup for the next layer
        currentLayerWidth = currentLayerWidth->next;
    }
    // +1 because depth is an index and thus one smaller than the
    // existing used size the array and we need one extra extra
    // element for the NaN terminator.
    CCE = realloc(CCE, sizeof(double) * (depth + 1));

    assert(CCE != NULL);

    // NaN terminate the array
    CCE[depth] = NAN;

    // ensure both queues are empty
    while(!isEmpty(currentLayer)) {
        dequeue(currentLayer);
    }
    while(!isEmpty(nextLayer)) {
        dequeue(nextLayer);
    }
    free(currentLayer);
    free(nextLayer);
    free(bfsLayer);

    return CCE;
}

void freellNode(void** llData){
    free((void*)*llData);
}

/*
    Fills memory pointed to by sequenceMemory with length random
    values between 0 and 4. Used for Testing.

    Input Args:
        int length: The length of the memory to fill.
    Output Args:
        int* sequenceMemory: The ponter to the memory to be filled.
*/
void genRandomSequence(int* sequenceMemory, int length) {
    int i;

    assert(sequenceMemory != NULL);

    for(i=0; i<length; ++i) {
        sequenceMemory[i] = random() % 5;
    }
}

int main(int argc, char const *argv[])
{
    int i;
    int* seq = malloc(sizeof(int) * 51);
    double* cces;

    srand(0);
    tRoot* root = createTree(5);
    for(i=0; i<10000; ++i) {
        genRandomSequence(seq, 50);
        insertSequence(root, seq, 50);
    }
    // for(i=0; i<50; ++i){
    //     insertSequence(root, seq+i, 50-i);
    // }
    for(i=0; i<50; ++i) {
    cces=calcCCEs(root);
    free(cces);
    }

    free(seq);
    
    freeTree(root);
    return 0;
}
