#ifndef HASHFILTER_H
#define HASHFILTER_H
#include <iostream>
#include <cstring>
#include <algorithm>

__always_inline uint32_t Get_Log(uint32_t x)
{
    uint32_t cnt = 0;
    for (; x > 0; cnt++)
        x >>= 1;
    return cnt;
}

__always_inline uint32_t Get_Size(uint32_t x)
{
    return x == 0 ? 0 : 1 << (Get_Log(x) - 1);
}

class HashFilter
{
private:
    bool *pool;
    uint32_t len;
    uint32_t radix;

public:
    HashFilter()
    {
    }
    ~HashFilter()
    {
        free(pool);
    }
    void Init(uint32_t n)
    {
        len = n;
        pool = (bool *)malloc(sizeof(bool) * len);
    }
    void Construct(uint32_t u)
    {
        constexpr int heuristic_factor = 32;
        uint32_t psize = std::max(Get_Size(u + 1), (uint32_t)4) * heuristic_factor;
        radix = psize - 1;
        pool = (bool *)malloc(sizeof(bool) * psize);
        memset(pool, false, sizeof(bool) * (radix + 1));
    }
    __always_inline void Set(uint32_t v)
    {
        pool[v & radix] = true;
    }
    __always_inline bool Try(uint32_t v)
    {
        return pool[v & radix];
    }
};

#endif