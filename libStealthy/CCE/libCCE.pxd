cdef extern from "linkedList.h":
    ctypedef struct llRoot:
        int length
        llNode* head
        llNode* tail
    ctypedef struct llNode:
        void* contents
        llNode* next
    llRoot* createLinkedList() nogil
    llNode* createllNode(void* value) nogil
    int isEmpty(llRoot* r) nogil
    void insertllBack(llRoot* r, void* value) nogil
    void insertllFront(llRoot* r, void* value) nogil
    void freeLL(llRoot* r) nogil
    void applyToList(llRoot* r, void (*func)(void**)) nogil
    void printList(llRoot* r, void (*printFunction)(void**)) nogil

cdef extern from "arenaAlloc.h":
    ctypedef struct arenaManager:
        pass
    void resetArena(arenaManager* m)

cdef extern from "ntree.h":
    ctypedef struct tNode:
        tNode** children
        unsigned long count
    ctypedef struct tRoot:
        tNode** children
        unsigned long branchingFactor;
        unsigned long layerCount;
        llRoot* layerWidth;
        arenaManager* arenaManagement;
    tNode* createNode(int branchingFactor)
    tRoot* createTree(int branchingFactor)
    void freeTree(tRoot* root) nogil
 
cdef extern from "CCE.h":
    void insertSequence(tRoot* root, int* sequence, int length)
    void printLayer(void** n)
    double* calcCCEs(tRoot* root)

cdef extern from "<float.h>":
    int isnan(double x) nogil