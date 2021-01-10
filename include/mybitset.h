#ifndef MYBITSET_H
#define MYBITSET_H
#include <iostream>
#include <cstring>
#include "tools.h"

class MyBitSet
{
private:
    uint64_t len;
    uint64_t *Buffer;

public:
    MyBitSet()
    {
    }
    void Init(uint32_t n)
    {
        len = n;
        Buffer = (uint64_t *)malloc(sizeof(uint64_t) * len);
        // memset(Buffer, 0, sizeof(uint64_t) * len);
    }
    ~MyBitSet()
    {
        free(Buffer);
    }
    __always_inline void Set(uint32_t Pos)
    {
        Buffer[Pos >> ShiftBits] |= (WordType)(1u) << (Pos & ModBits);
    }
    __always_inline bool Get(uint32_t Pos)
    {
        return Buffer[Pos >> ShiftBits] & ((WordType)(1u) << (Pos & ModBits));
    }
    __always_inline void SetWord(uint32_t Pos, WordType Word)
    {
        Buffer[Pos] = Word;
    }
    __always_inline
        WordType
        GetWord(uint32_t Pos)
    {
        return Buffer[Pos];
    }
    __always_inline void Clr()
    {
        memset(Buffer, 0, sizeof(uint64_t) * len);
    }
};

#endif