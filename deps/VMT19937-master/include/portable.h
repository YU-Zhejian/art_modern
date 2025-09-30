#pragma once

// #include "simd_config.h"

#include <cstdint>
#include <cstddef>

#ifdef _MSC_VER

inline int64_t popcnt(uint64_t x)
{
    return __popcnt64(x);
}

inline int32_t popcnt(uint32_t x)
{
    return __popcnt(x);
}

#else

inline int64_t popcnt(uint64_t x)
{
    return __builtin_popcountll(x);
}

inline int32_t popcnt(uint32_t x)
{
    return __builtin_popcountl(x);
}

#endif
