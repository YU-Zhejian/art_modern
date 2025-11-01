

#include "libam_support/seq/str_to_qual.h"

#include "libam_support/Constants.h"
#include "libam_support/Dtypes.h"

// NOLINTBEGIN
#if defined(__SSE2__) || defined(__AVX2__) || defined(__MMX__)
#include <immintrin.h> // FIXME: Change to <x86intrin.h> or <intrin.h> (MSVC) to improve portability
#endif
// NOLINTEND

#include <stdlib.h>

void str_to_qual_for_loop(am_qual_t* qual, const char* str, const size_t qlen)
{
    for (size_t i = 0; i < qlen; ++i) {
        qual[i] = (char)(str[i] - AM_PHRED_OFFSET);
    }
}
void str_to_qual_avx2(am_qual_t* qual, const char* str, const size_t qlen)
{
    size_t i = 0;
#ifdef __AVX2__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 32; // AVX2 processes 32 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 5) << 5; // Align to 32-byte boundary
    __m256i phred_offset_vec = _mm256_set1_epi8((uint8_t)(AM_PHRED_OFFSET));
    for (; i < aligned_size; i += num_elements_per_simd) {
        __m256i str_vec = _mm256_loadu_si256((const __m256i*)(&str[i]));
        __m256i result_vec = _mm256_sub_epi8(str_vec, phred_offset_vec);
        _mm256_storeu_si256((__m256i*)(&qual[i]), result_vec);
    }
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        qual[i] = (char)(str[i] - AM_PHRED_OFFSET);
    }
}
void str_to_qual_sse2(am_qual_t* qual, const char* str, const size_t qlen)
{
    size_t i = 0;
#ifdef __SSE2__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 16; // SSE2 processes 16
    const size_t aligned_size = (qlen >> 4) << 4; // Align to 16-byte boundary
    __m128i phred_offset_vec = _mm_set1_epi8((uint8_t)(AM_PHRED_OFFSET));
    for (; i < aligned_size; i += num_elements_per_simd) {
        __m128i str_vec = _mm_loadu_si128((const __m128i*)(&str[i]));
        __m128i result_vec = _mm_sub_epi8(str_vec, phred_offset_vec);
        _mm_storeu_si128((__m128i*)(&qual[i]), result_vec);
    }
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        qual[i] = (char)(str[i] - AM_PHRED_OFFSET);
    }
}
void str_to_qual_mmx(am_qual_t* qual, const char* str, const size_t qlen)
{
    size_t i = 0;
#ifdef __MMX__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 8; // MMX processes 8 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 3) << 3; // Align to 8-byte boundary
    __m64 phred_offset_vec = _mm_set1_pi8((uint8_t)(AM_PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m64 str_vec = *(const __m64*)(&str[i]);
        __m64 result_vec = _mm_sub_pi8(str_vec, phred_offset_vec);
        *(__m64*)(&qual[i]) = result_vec;
    }
    _mm_empty(); // Empty the MMX state
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        qual[i] = (char)(str[i] - AM_PHRED_OFFSET);
    }
}
void str_to_qual_comb(am_qual_t* qual, const char* str, size_t qlen)
{
    size_t i = 0;
#ifdef __AVX2__
    // Process with AVX2 (32 elements at a time)
    const size_t avx2_aligned_size = (qlen >> 5) << 5; // Align to 32-byte boundary
    if (avx2_aligned_size > 0) {
        const __m256i phred_offset_vec = _mm256_set1_epi8((uint8_t)(AM_PHRED_OFFSET));
        for (; i < avx2_aligned_size; i += 32) {
            const __m256i str_vec = _mm256_loadu_si256((const __m256i*)(&str[i]));
            const __m256i result_vec = _mm256_sub_epi8(str_vec, phred_offset_vec);
            _mm256_storeu_si256((__m256i*)(&qual[i]), result_vec);
        }
    }
#endif
#ifdef __SSE2__
    // Process remaining with SSE2 (16 elements at a time)
    const size_t sse2_aligned_size = (qlen >> 4) << 4; // Align to 16-byte boundary
    if (i < sse2_aligned_size) {
        const __m128i phred_offset_vec = _mm_set1_epi8((uint8_t)(AM_PHRED_OFFSET));
        for (; i < sse2_aligned_size; i += 16) {
            const __m128i str_vec = _mm_loadu_si128((const __m128i*)(&str[i]));
            const __m128i result_vec = _mm_sub_epi8(str_vec, phred_offset_vec);
            _mm_storeu_si128((__m128i*)(&qual[i]), result_vec);
        }
    }
#endif
#ifdef __MMX__
    // Process remaining with MMX (8 elements at a time)
    const size_t mmx_aligned_size = (qlen >> 3) << 3; // Align to 8-byte boundary
    if (i < mmx_aligned_size) {
        const __m64 phred_offset_vec = _mm_set1_pi8((uint8_t)(AM_PHRED_OFFSET));
        for (; i < mmx_aligned_size; i += 8) {
            const __m64 str_vec = *(const __m64*)(&str[i]);
            const __m64 result_vec = _mm_sub_pi8(str_vec, phred_offset_vec);
            *(__m64*)(&qual[i]) = result_vec;
        }
        _mm_empty(); // Empty the MMX state
    }
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        qual[i] = (char)(str[i] - AM_PHRED_OFFSET);
    }
}
