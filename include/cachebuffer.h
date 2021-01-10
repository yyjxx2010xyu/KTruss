#ifndef CACHEBUFFER_H
#define CACHEBUFFER_H
#include <iostream>
#include <cstring>

//  Reference https://blog.cheyulin.me/?p=840
class CacheBuffer
{
    uint32_t *Local, *Origin;
    uint32_t bufferSize, bufferLen;
    volatile uint32_t *originPtr;

    void pushBuffer()
    {
        uint32_t ptr = __sync_fetch_and_add(originPtr, bufferSize);
        memcpy(Origin + ptr, Local, sizeof(uint32_t) * bufferSize);
        bufferSize = 0;
    }

public:
    CacheBuffer(){};
    void Init(uint32_t _bufferLen)
    {
        bufferLen = _bufferLen;
        Local = (uint32_t *)malloc(sizeof(uint32_t) * (_bufferLen + 1));
    }
    void Set(uint32_t *_Origin, volatile uint32_t *_originPtr)
    {
        Origin = _Origin;
        originPtr = _originPtr;
        bufferSize = 0;
    }
    ~CacheBuffer()
    {
        free(Local);
    }
    void push(uint32_t data)
    {
        Local[bufferSize++] = data;
        if (bufferSize >= bufferLen)
            pushBuffer();
    }

    void pushBufferAnyway()
    {
        if (bufferSize != 0)
            pushBuffer();
    }
};

#endif