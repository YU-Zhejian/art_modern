/**
 * @file      qreverse.h
 * @brief     Functions to reverse sequences of various element sizes, re-implemented in C.
 * @copyright
 * 			2017 Wunkolo
 * 			2025 YU Zhejian
 * @license   MIT
 */

#include "art_modern_config.h" // NOLINT: for AM_NO_Q_REVERSE

#if !defined(AM_NO_Q_REVERSE)

#ifndef ART_MODERN_LIBAM_SUPPORT_SEQ_QREVERSE_H
#define ART_MODERN_LIBAM_SUPPORT_SEQ_QREVERSE_H

#include "libam_support/utils/c_macro_utils.h"

ART_MODERN_BEGIN_C_DECLS

#include <stddef.h>
#include <stdint.h>

// x86
#if defined(__i386__) || defined(__x86_64__)
#if defined(_MSC_VER)
#include <intrin.h>
#define Swap64 _byteswap_uint64
#define Swap32 _byteswap_ulong
#define Swap16 _byteswap_ushort

#elif defined(__GNUC__) || defined(__clang__)
#include <x86intrin.h>

#define Swap64 __builtin_bswap64
#define Swap32 __builtin_bswap32
#define Swap16 __builtin_bswap16
#endif

// ARM
#elif defined(__ARM_NEON)
#include <arm_neon.h>
#if defined(_MSC_VER)
#define Swap64 _byteswap_uint64
#define Swap32 _byteswap_ulong
#define Swap16 _byteswap_ushort

#elif defined(__GNUC__) || defined(__clang__)
#define Swap64 __builtin_bswap64
#define Swap32 __builtin_bswap32
#define Swap16 __builtin_bswap16

#endif

// Pure
#else

inline uint64_t Swap64(uint64_t x)
{
    return (((x & 0x00000000000000FF) << 56) | ((x & 0x000000000000FF00) << 40) | ((x & 0x0000000000FF0000) << 24)
        | ((x & 0x00000000FF000000) << 8) | ((x & 0x000000FF00000000) >> 8) | ((x & 0x0000FF0000000000) >> 24)
        | ((x & 0x00FF000000000000) >> 40) | ((x & 0xFF00000000000000) >> 56));
}

inline uint32_t Swap32(uint32_t x)
{
    return (((x & 0x000000FF) << 24) | ((x & 0x0000FF00) << 8) | ((x & 0x00FF0000) >> 8) | ((x & 0xFF000000) >> 24));
}

inline uint16_t Swap16(uint16_t x) { return (((x & 0x00FF) << 8) | ((x & 0xFF00) >> 8)); }

#endif

// One byte elements
void qReverse_1(void* Array, size_t Count);

// Two byte elements
void qReverse_2(void* Array, size_t Count);

// Four byte elements
void qReverse_4(void* Array, size_t Count);

// 8 byte elements
void qReverse_8(void* Array, size_t Count);

ART_MODERN_END_C_DECLS

#endif
#endif
