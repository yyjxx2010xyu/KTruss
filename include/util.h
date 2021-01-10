#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <omp.h>

class uvList
{
public:
    uint32_t u;
    uint32_t v;
    uint32_t uvId;
    uvList (uint32_t _u, uint32_t _v, uint32_t _uvId) : u(_u), v(_v), uvId(_uvId) {} 
};

bool orderComp (const uint32_t& u, const uint32_t& v);


struct Edge {
    uint32_t u;
    uint32_t v;
} __attribute__ ((aligned (4)));

bool operator <(const Edge &e1, const Edge &e2);

void parmemset(void *const ptr, int const c, uint32_t const len);


#endif