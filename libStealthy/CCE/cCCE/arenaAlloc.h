#pragma once
#include <stdint.h>

typedef struct arenaStruct {
    uint8_t* arena;
    unsigned long currentOffset;
} arena;

typedef struct arenaMetadata {
    arena* arenas;
    unsigned long arenaSize;
    int currentArena;
    int arenaCount;
} arenaManager;


/*
    Returns a pointer to an arenaManager structure which holds meta-
    data reguarding the size, number and in memory positions of the
    allocation arenas.  It also reserves virtual memory for the first
    arena equal to arenaSize.
    
    Input Args:
        unsigned long arenaSize: The size of the pre-allocation arenas
            once one is full another will be allocated
    Returns:
        arenaManager*: A pointer to an arena manager 
*/
arenaManager* setupArena(unsigned long arenaSize);

/*
    Returns a pointer to memory inside an arena handled by the
    arenaManager m.  If the size of the requested memory is larger
    than what the arena has left a new arena will be created and any
    memory left in the orignal arena will not be used.

    Any memory returned by this is not protected by OS level
    segmentation and as a result care should be taken to ensure one
    doesn't trample on the memory of another call.

    Input Args:
        arenaManager* m: A pointer to an arenaManager.
        unsigned long size: the size in bytes of the memory you would
            like set set aside for you.
    Returns:
        void*: A pointer to memory in the arena.
*/
void* arenaAlloc(arenaManager* m, unsigned long size);

/*
    Frees all memory assocated wih an arenaManager and then frees the
    arenaManager itself.
    
    Input Args:
        arenaManager* m: A pointer to an arenaManager.
*/
void teardownArena(arenaManager* m);

/*
    Resets the currentArena index allowing for re-use of pre-
    allocated memory.  In certan use cases this can provide a massive
    speedup over simply freeing all arenas and re-alocating them.

    Input Args:
        arenaManager* m: A pointer to an arenaManager.
*/
void resetArena(arenaManager* m);