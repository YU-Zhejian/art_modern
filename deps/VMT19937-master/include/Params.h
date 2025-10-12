#pragma once

#include <cstdint>
#include <cstddef>

namespace Details {

struct MT19937Params
{
    static const size_t s_nBits = 19937;
    static const size_t s_wordSizeBits = 32;
    static const int s_N = s_nBits / s_wordSizeBits + (s_nBits % s_wordSizeBits != 0);  // 624
    static const int s_M = 397;
    const static size_t s_nMatrixBits = s_nBits;

    static const uint32_t s_temperMask1 = 0x9d2c5680UL;
    static const uint32_t s_temperMask2 = 0xefc60000UL;

    static const uint32_t s_matrixA = 0x9908b0dfUL;   // constant vector a
    static const uint32_t s_upperMask = 0x80000000UL; // most significant w-r bits
    static const uint32_t s_lowerMask = 0x7fffffffUL; // least significant r bits
};

struct SFMT19937Params
{
    static const size_t s_nBits = 19937;
    
    static const size_t s_wordSizeBits = 128;
    static const int s_N = s_nBits / s_wordSizeBits + (s_nBits % s_wordSizeBits != 0);  // 156
    static const int s_M = 122;
    static const size_t s_nMatrixBits = s_N * s_wordSizeBits;

    static const uint32_t s_SFMT_MSK1 = 0xdfffffefU;
    static const uint32_t s_SFMT_MSK2 = 0xddfecb7fU;
    static const uint32_t s_SFMT_MSK3 = 0xbffaffffU;
    static const uint32_t s_SFMT_MSK4 = 0xbffffff6U;
};

};
