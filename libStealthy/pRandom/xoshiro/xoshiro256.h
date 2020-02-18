/*
* U.S. Naval Research Laboratory
* Center for High Assurance Computer Systems
*/

#ifndef XOSHIRO256
#define XOSHIRO256 1
#include <stdint.h>
static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}
uint64_t next(uint64_t* s);
void jump(uint64_t* s);
void long_jump(uint64_t* s);
#endif