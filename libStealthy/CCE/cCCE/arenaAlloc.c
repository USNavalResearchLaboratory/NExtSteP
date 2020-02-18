#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "arenaAlloc.h"

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
arenaManager* setupArena(unsigned long arenaSize) {
    arenaManager* m;

    assert(arenaSize > 0);

    m = (arenaManager*)malloc(sizeof(arenaManager));
    assert(m != NULL);

    m->arenas = (arena*)malloc(sizeof(arena));
    assert(m->arenas != NULL);

    m->arenas[0].arena = (uint8_t*)malloc(arenaSize);
    assert(m->arenas[0].arena != NULL);

    m->arenas[0].currentOffset = 0;
    m->arenaCount = 1;
    m->arenaSize = arenaSize;
    m->currentArena = 0;

    return m;
}

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
void* arenaAlloc(arenaManager* m, unsigned long size) {
    arena* currentArena;
    unsigned long previousOffset;

    assert(m != NULL);
    assert(size <= m->arenaSize);
    assert(m->arenas != NULL);

    currentArena = &(m->arenas[m->currentArena]);

    
    if((currentArena->currentOffset + size) >= m->arenaSize) {
        m->currentArena += 1;
        if(m->currentArena >= m->arenaCount) {
            #ifdef DEBUG_ALLOC
                printf("%lu %lu\n", (m->arenas[m->arenaCount-1].currentOffset + size), m->arenaSize);
            #endif
            m->arenaCount += 1;
            m->arenas = realloc(m->arenas, m->arenaCount * sizeof(arena));
            assert(m->arenas != NULL);

            currentArena = &(m->arenas[m->currentArena]);
            currentArena->arena = (uint8_t*)malloc(m->arenaSize);
            assert(currentArena->arena != NULL);

        } else {
            currentArena = &(m->arenas[m->currentArena]);
        }
        currentArena->currentOffset = 0;

    }
    previousOffset = currentArena->currentOffset;
    currentArena->currentOffset += size;

    // Get the ponter to the start of the space requested.
    return (void*)(&(currentArena->arena[previousOffset]));
}

/*
    Frees all memory assocated wih an arenaManager and then frees the
    arenaManager itself.

    Input Args:
        arenaManager* m: A pointer to an arenaManager.
*/
void teardownArena(arenaManager* m) {
    int i;
    assert(m != NULL);
    assert(m->arenas != NULL);

    for(i = 0; i < m->arenaCount; ++i) {
        assert(m->arenas[i].arena != NULL);
        free(m->arenas[i].arena);
    }
    free(m->arenas);
    free(m);
}

/*
    Resets the currentArena index allowing for re-use of pre-
    allocated memory.  In certan use cases this can provide a massive
    speedup over simply freeing all arenas and re-alocating them.

    Input Args:
        arenaManager* m: A pointer to an arenaManager.
*/
void resetArena(arenaManager* m) {
    m->currentArena = 0;
}

// int main(int argc, char const *argv[])
// {
//     arenaManager* m = setupArena(5);
//     uint16_t* p = arenaAlloc(m, sizeof(uint16_t));
//     uint16_t* d = arenaAlloc(m, sizeof(uint16_t));
//     uint16_t* e = arenaAlloc(m, sizeof(uint16_t));
//     printf("reset\n");
//     resetArena(m);
//     p = arenaAlloc(m, sizeof(uint16_t));
//     d = arenaAlloc(m, sizeof(uint16_t));
//     e = arenaAlloc(m, sizeof(uint16_t));
//     p = arenaAlloc(m, sizeof(uint16_t));
//     d = arenaAlloc(m, sizeof(uint16_t));
//     e = arenaAlloc(m, sizeof(uint16_t));
//     printf("reset\n");
//     resetArena(m);
//     p = arenaAlloc(m, sizeof(uint16_t));
//     d = arenaAlloc(m, sizeof(uint16_t));
//     e = arenaAlloc(m, sizeof(uint16_t));
//     p = arenaAlloc(m, sizeof(uint16_t));
//     d = arenaAlloc(m, sizeof(uint16_t));
//     e = arenaAlloc(m, sizeof(uint16_t));
//     p = arenaAlloc(m, sizeof(uint16_t));
//     d = arenaAlloc(m, sizeof(uint16_t));
//     e = arenaAlloc(m, sizeof(uint16_t));
//     p = arenaAlloc(m, sizeof(uint16_t));
//     d = arenaAlloc(m, sizeof(uint16_t));
//     e = arenaAlloc(m, sizeof(uint16_t));
//     printf("end\n");

//     printf("%p\n", p);
//     *p = (uint16_t)5;
//     printf("%u\n", *p);
//     teardownArena(m);
//     return 0;
// }