#ifndef LIBPOPC_H
#define LIBPOPCNT_H
#if 0
//  http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
uint64_t mypopcount(uint64_t x)
{
  const static uint64_t m1 = 0x5555555555555555ll;
  const static uint64_t m2 = 0x3333333333333333ll;
  const static uint64_t m4 = 0x0F0F0F0F0F0F0F0Fll;
  const static uint64_t h01 = 0x0101010101010101ll;

  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}
#endif

//  内嵌汇编，似乎会快一点！ 11.17 但是需要SSE4.2支持！
__always_inline
    uint64_t
    inspopcount(const uint64_t x)
{
  uint64_t cnt = 0;
  __asm__ __volatile__(
      "popcnt %1, %1  \n\t"
      "mov %1, %0     \n\t"
      : "+r"(cnt)
      : "r"(x));
  return cnt;
}

#endif
