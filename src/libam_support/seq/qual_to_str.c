/**
 * Copyright 2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/
#include "libam_support/seq/qual_to_str.h"
#include "libam_support/utils/c_macro_utils.h"

#include "libam_support/Constants.h"
#include "libam_support/Dtypes.h"

#ifdef INCLUDE_INTEL_SIMD
#include INCLUDE_INTEL_SIMD // NOLINT
#endif

#include <stdlib.h>

void qual_to_str_mmx(const am_qual_t* qual, char* str, const size_t qlen)
{
    size_t i = 0;
#ifdef __MMX__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 8; // MMX processes 8 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 3) << 3; // Align to 8-byte boundary
    __m64 phred_offset_vec = _mm_set1_pi8((uint8_t)(AM_PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m64 qual_vec = *(const __m64*)(&qual[i]);
        __m64 result_vec = _mm_add_pi8(qual_vec, phred_offset_vec);
        *(__m64*)(&str[i]) = result_vec;
    }
    _mm_empty(); // Empty the MMX state
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        str[i] = (char)(qual[i] + AM_PHRED_OFFSET);
    }
}

void qual_to_str_sse2(const am_qual_t* qual, char* str, const size_t qlen)
{
    size_t i = 0;
#ifdef __SSE2__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 16; // SSE2 processes 16 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 4) << 4; // Align to 16-byte boundary
    __m128i phred_offset_vec = _mm_set1_epi8((uint8_t)(AM_PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m128i qual_vec = _mm_loadu_si128((const __m128i*)(&qual[i]));
        __m128i result_vec = _mm_add_epi8(qual_vec, phred_offset_vec);
        _mm_storeu_si128((__m128i*)(&str[i]), result_vec);
    }
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        str[i] = (char)(qual[i] + AM_PHRED_OFFSET);
    }
}

void qual_to_str_avx2(const am_qual_t* qual, char* str, const size_t qlen)
{
    size_t i = 0;
#ifdef __AVX2__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 32; // AVX2 processes 32 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 5) << 5; // Align to 32-byte boundary
    __m256i phred_offset_vec = _mm256_set1_epi8((uint8_t)(AM_PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m256i qual_vec = _mm256_loadu_si256((const __m256i*)(&qual[i]));
        __m256i result_vec = _mm256_add_epi8(qual_vec, phred_offset_vec);
        _mm256_storeu_si256((__m256i*)(&str[i]), result_vec);
    }
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        str[i] = (char)(qual[i] + AM_PHRED_OFFSET);
    }
}

void qual_to_str_comb(const am_qual_t* qual, char* str, const size_t qlen)
{
    size_t i = 0;

#ifdef __AVX2__
    // Process with AVX2 (32 elements at a time)
    const size_t avx2_aligned_size = (qlen >> 5) << 5; // Align to 32-byte boundary
    if (avx2_aligned_size > 0) {
        const __m256i phred_offset_vec = _mm256_set1_epi8((uint8_t)(AM_PHRED_OFFSET));
        for (; i < avx2_aligned_size; i += 32) {
            const __m256i qual_vec = _mm256_loadu_si256((const __m256i*)(&qual[i]));
            const __m256i result_vec = _mm256_add_epi8(qual_vec, phred_offset_vec);
            _mm256_storeu_si256((__m256i*)(&str[i]), result_vec);
        }
    }
#endif

#ifdef __SSE2__
    // Process remaining with SSE2 (16 elements at a time)
    const size_t sse2_aligned_size = (qlen >> 4) << 4; // Align to 16-byte boundary
    if (i < sse2_aligned_size) {
        const __m128i phred_offset_vec = _mm_set1_epi8((uint8_t)(AM_PHRED_OFFSET));
        for (; i < sse2_aligned_size; i += 16) {
            const __m128i qual_vec = _mm_loadu_si128((const __m128i*)(&qual[i]));
            const __m128i result_vec = _mm_add_epi8(qual_vec, phred_offset_vec);
            _mm_storeu_si128((__m128i*)(&str[i]), result_vec);
        }
    }
#endif

#ifdef __MMX__
    // Process remaining with MMX (8 elements at a time)
    const size_t mmx_aligned_size = (qlen >> 3) << 3; // Align to 8-byte boundary
    if (i < mmx_aligned_size) {
        const __m64 phred_offset_vec = _mm_set1_pi8((uint8_t)(AM_PHRED_OFFSET));
        for (; i < mmx_aligned_size; i += 8) {
            const __m64 qual_vec = *(const __m64*)(&qual[i]);
            const __m64 result_vec = _mm_add_pi8(qual_vec, phred_offset_vec);
            *(__m64*)(&str[i]) = result_vec;
        }
        _mm_empty(); // Empty the MMX state
    }
#endif

    // Handle the remaining elements that do not fit into any SIMD register
    for (; i < qlen; ++i) {
        str[i] = (char)(qual[i] + AM_PHRED_OFFSET);
    }
}

void qual_to_str_for_loop(const am_qual_t* qual, char* str, const size_t qlen)
{
    // Handle the remaining elements that do not fit into a full SIMD register
    for (size_t i = 0; i < qlen; ++i) {
        str[i] = (char)(qual[i] + AM_PHRED_OFFSET);
    }
}
