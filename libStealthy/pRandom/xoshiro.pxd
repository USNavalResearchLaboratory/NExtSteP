from libc.stdint cimport uint64_t

cdef extern from "xoshiro256.h":
    uint64_t next(uint64_t* s)
    void jump(uint64_t* s)
    void long_jump(uint64_t* s)