#pragma once

#include <cstdint>
#include <cstddef>

#include "macros.h"
#if defined(_MSC_VER) && (_M_IX86_FP==2 || defined(_M_X64))
#  define __SSE2__
#  define __SSE4_1__
#  define __SSE4_2__
#endif

#if defined(__AVX512F__)
#   define SIMD_N_BITS 512
#elif defined(__AVX2__)
#   define SIMD_N_BITS 256
#elif defined(__AVX__)
#   error AVX2 is needed
#elif defined(__SSE4_2__)
#   define SIMD_N_BITS 128
#elif defined(__SSE2__)
#   error SSE4.2 is needed
#endif

#ifdef SIMD_N_BITS
#   include <immintrin.h>
#   ifdef _MSC_VER
#       include <intrin.h>
#   endif
#else
#error "SIMD_N_BITS not defined"
#endif

