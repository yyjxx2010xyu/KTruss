#ifndef FIXEDVECTOR_H
#define FIXEDVECTOR_H

#include <iostream>

// Just For Competition, Fixed Size Vector, Provide atomic push_back
class FixedVector
{
private:
    uint32_t _len;
    uint32_t _size;
    uint32_t * arr;
public:
    FixedVector(uint32_t len)
    {
        _len = len;
        _size = 0;
        arr = (uint32_t *)malloc(sizeof(uint32_t) * _len);
    }
    
    ~FixedVector()
    {
        free(arr);
    }
    
    uint32_t size()
    {
        return _size;
    }
    void push_back(uint32_t v)
    {
        uint32_t ptr;
        do
        {
            ptr = _size;
        } while (!__sync_bool_compare_and_swap(&_size, ptr, ptr + 1));
        arr[ptr] = v;
    }
    uint32_t operator [](const uint32_t p)
	{
		return arr[p];
	}
    void clear()
    {
        _size = 0;
    }
    void swap(FixedVector & Other)
    {
        std::swap(this->_size, Other._size);
        std::swap(this->arr, Other.arr);
    }

    uint32_t * GetSizePtr()
    {
        return &_size;
    }
    uint32_t * GetArrPtr()
    {
        return arr;
    }
};
#endif