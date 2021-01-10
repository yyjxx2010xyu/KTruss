#ifndef TOOLS_H
#define TOOLS_H
#include <omp.h>
// #define TIMER_ON
// #define TRI_ON
#define REV_ON
// #define FILTER_ON
// #define DEBUG_ON
// #define ASSERT_ON
#include <iostream>
constexpr uint32_t ENDLIST = 0x3f3f3f3f;
using WordType = uint64_t;
constexpr uint32_t WordBits = sizeof(WordType) * 8;
constexpr uint32_t ShiftBits = 6;
constexpr uint32_t ModBits = 0x0000003f;
constexpr uint32_t MINUS = 0x7fffffff;
constexpr uint32_t EDGEDELETE = 0x7fffffff;
constexpr uint32_t INF = 0x7fffffff;
constexpr uint32_t CLRT = 1000;
constexpr uint32_t BUFFER_SIZE = 1 << 12;
static uint32_t ThreadNum  = omp_get_num_procs();

#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

#if !defined(__likely)
#define __likely(...) (__builtin_expect((__VA_ARGS__), 1))
#endif // !defined(__likely)

#if !defined(__unlikely)
#define __unlikely(...) (__builtin_expect((__VA_ARGS__), 0))
#endif // !defined(__unlikely)


#endif