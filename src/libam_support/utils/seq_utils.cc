/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "seq_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"

// NOLINTBEGIN
#if defined(__SSE2__) || defined(__AVX2__) || defined(__MMX__)
#include <immintrin.h>
#endif
// NOLINTEND

#include <htslib/sam.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace labw::art_modern {
constexpr char rev_comp_trans_2[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
    23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73, 74, 75, 76, 77, 78, 79, 80,
    81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 };
/**
 *
 * Generated using:
 *
 * ```python
 * ",".join([str(i) if chr(i) in 'AGCT' else str(ord(chr(i).upper())) if chr(i) in 'agct' else str(ord('N'))for i in
 * range(0, 128)])
 * ```
 */
constexpr char normalization_matrix[] = { 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 65, 78, 67, 78, 78, 78, 71, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 84, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 65, 78, 67, 78, 78, 78, 71, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 84, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78 };

std::string qual_to_str_mmx(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    size_t i = 0;
#ifdef __MMX__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 8; // MMX processes 8 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 3) << 3; // Align to 8-byte boundary
    __m64 phred_offset_vec = _mm_set1_pi8(static_cast<uint8_t>(PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m64 qual_vec = *reinterpret_cast<const __m64*>(&qual[i]);
        __m64 result_vec = _mm_add_pi8(qual_vec, phred_offset_vec);
        *reinterpret_cast<__m64*>(&retq[i]) = result_vec;
    }
    _mm_empty(); // Empty the MMX state
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        retq[i] = static_cast<char>(qual[i] + PHRED_OFFSET);
    }
    return retq;
}

std::string qual_to_str_sse2(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    size_t i = 0;
#ifdef __SSE2__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 16; // SSE2 processes 16 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 4) << 4; // Align to 16-byte boundary
    __m128i phred_offset_vec = _mm_set1_epi8(static_cast<uint8_t>(PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m128i qual_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&qual[i]));
        __m128i result_vec = _mm_add_epi8(qual_vec, phred_offset_vec);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(&retq[i]), result_vec);
    }
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        retq[i] = static_cast<char>(qual[i] + PHRED_OFFSET);
    }
    return retq;
}

std::string qual_to_str_avx2(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    size_t i = 0;
#ifdef __AVX2__
    // NOLINTBEGIN
    const size_t num_elements_per_simd = 32; // AVX2 processes 32 uint8_t elements at a time
    const size_t aligned_size = (qlen >> 5) << 5; // Align to 32-byte boundary
    __m256i phred_offset_vec = _mm256_set1_epi8(static_cast<uint8_t>(PHRED_OFFSET));

    for (; i < aligned_size; i += num_elements_per_simd) {
        __m256i qual_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&qual[i]));
        __m256i result_vec = _mm256_add_epi8(qual_vec, phred_offset_vec);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(&retq[i]), result_vec);
    }
    // NOLINTEND
#endif
    // Handle the remaining elements that do not fit into a full SIMD register
    for (; i < qlen; ++i) {
        retq[i] = static_cast<char>(qual[i] + PHRED_OFFSET);
    }
    return retq;
}

std::string qual_to_str_for_loop(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    // Handle the remaining elements that do not fit into a full SIMD register
    for (size_t i = 0; i < qlen; ++i) {
        retq[i] = static_cast<char>(qual[i] + PHRED_OFFSET);
    }
    return retq;
}

std::string qual_to_str_foreach(const am_qual_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    std::memcpy(retq.data(), qual, qlen);
    std::for_each(retq.begin(), retq.end(), [](char& c) { c += PHRED_OFFSET; });
    return retq;
}

std::string qual_to_str(const am_qual_t* qual, const size_t qlen)
{
#if defined(__MMX__) || defined(__SSE2__) || defined(__AVX2__)
    if (qlen <= 100) {
        return qual_to_str_mmx(qual, qlen);
    }
    if (qlen <= 200) {
        return qual_to_str_sse2(qual, qlen);
    }
    return qual_to_str_avx2(qual, qlen);
#else
    return qual_to_str_foreach(qual, qlen);
#endif
}

std::string qual_to_str(const std::vector<am_qual_t>& qual) { return qual_to_str(qual.data(), qual.size()); }

std::string comp(const std::string& dna)
{
    std::string rets;
    rets.resize(dna.length());
    for (decltype(dna.length()) i = 0; i < dna.length(); i++) {
        rets[i] = rev_comp_trans_2[dna[i] & 0xFF];
    }
    return rets;
}
std::string revcomp(const std::string& dna)
{
    auto rets = comp(dna);
    reverse(rets.begin(), rets.end());
    return rets;
}

[[maybe_unused]] std::string normalize(const std::string& dna)
{
    std::string rets = dna;
    std::for_each(rets.begin(), rets.end(), [](char& c) { c = normalization_matrix[c & 0xFF]; });
    return rets;
}

std::string cigar_arr_to_str(const std::vector<am_cigar_t>& cigar_arr)
{
    return cigar_arr_to_str_optim(cigar_arr.data(), cigar_arr.size());
}

std::string cigar_arr_to_str_old(const am_cigar_t* cigar_arr, const size_t n)
{
    std::ostringstream oss;
    for (size_t i = 0; i < n; i += 1) {
        oss << (cigar_arr[i] >> BAM_CIGAR_SHIFT);
        oss << BAM_CIGAR_STR[cigar_arr[i] & BAM_CIGAR_MASK];
    }
    return oss.str();
}

std::string cigar_arr_to_str_optim(const am_cigar_t* cigar_arr, const size_t n)
{
    std::string result;
    result.reserve(n << 1); // Rough estimate: each CIGAR op is at least 2 chars
    for (size_t i = 0; i < n; ++i) {
        result += std::to_string(cigar_arr[i] >> BAM_CIGAR_SHIFT);
        result += BAM_CIGAR_STR[cigar_arr[i] & BAM_CIGAR_MASK];
    }
    result.shrink_to_fit();
    return result;
}

void comp_inplace(std::string& dna)
{
    std::for_each(dna.begin(), dna.end(), [](char& c) { c = rev_comp_trans_2[c & 0xFF]; });
}

void revcomp_inplace(std::string& dna)
{
    comp_inplace(dna);
    reverse(dna.data(), dna.size());
}

void normalize_inplace(std::string& dna)
{
    std::for_each(dna.begin(), dna.end(), [](char& c) { c = normalization_matrix[c & 0xFF]; });
}

bool ends_with(const std::string& str, const std::string& suffix)
{
    if (str.length() < suffix.length()) {
        return false;
    }
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

} // namespace labw::art_modern