from libc.stdint cimport uint32_t, uint64_t
cdef extern from "pcg_basic.h":
    cdef struct pcg_state_setseq_64:
        uint64_t state
        uint64_t inc
    ctypedef pcg_state_setseq_64 pcg32_random_t
    void pcg32_srandom(uint64_t initstate, uint64_t initseq)
    void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
    uint32_t pcg32_random()
    uint32_t pcg32_random_r(pcg32_random_t* rng)
    uint32_t pcg32_boundedrand(uint32_t bound)
    uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)