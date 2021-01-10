#include <iostream>
#include "util.h"
#include <cstring>

bool operator<(const Edge &e1, const Edge &e2)
{
    return e1.u < e2.v || (e1.u == e2.u && e1.v < e2.v);
}

void parmemset(void *const ptr, int const c, uint32_t const len)
{
#pragma omp parallel
    {
        const uint32_t nthreads = omp_get_num_threads();
        const uint32_t tid = omp_get_thread_num();

        uint32_t const curlen = (len + nthreads - 1) / nthreads;
        uint32_t const beg = std::min(curlen * tid, len);
        uint32_t const end = std::min(beg + curlen, len);

        char *const cptr = (char *)ptr;
        memset(cptr + beg, c, end - beg);
    }

}