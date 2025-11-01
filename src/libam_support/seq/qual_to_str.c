#include "libam_support/seq/qual_to_str.h"

#include "libam_support/Constants.h"
#include "libam_support/Dtypes.h"

// NOLINTBEGIN
#if defined(__SSE2__) || defined(__AVX2__) || defined(__MMX__)
#include <immintrin.h>
#endif
// NOLINTEND

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


void qual_to_str_for_loop(const am_qual_t* qual, char* str, const size_t qlen)
{
    // Handle the remaining elements that do not fit into a full SIMD register
    for (size_t i = 0; i < qlen; ++i) {
        str[i] = (char)(qual[i] + AM_PHRED_OFFSET);
    }
}
