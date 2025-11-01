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

#include "libam_support/utils/seq_utils.hh"

#include "libam_support/seq/qual_to_str.h"
#include "libam_support/seq/str_to_qual.h"

#include "libam_support/Dtypes.h"

#include "ceu_check/ceu_check_cc_macro.h"

// Some compiler's -O3 optimization is better than ours
#if (defined(CEU_COMPILER_IS_INTEL_CLANG) || defined(CEU_COMPILER_IS_ICC) || defined(CEU_COMPILER_IS_CLANG)            \
    || defined(CEU_COMPILER_IS_GCC))                                                                                   \
    && !(defined(CEU_CM_IS_DEBUG))
#define IN_COMPILER_WE_TRUST
#endif

#include <htslib/sam.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
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

void str_to_qual(std::vector<am_qual_t>& qual, const std::string& str)
{
    str_to_qual(qual.data(), str.c_str(), str.size());
}

void str_to_qual(am_qual_t* qual, const char* str, size_t qlen)
{
#if !defined(IN_COMPILER_WE_TRUST)
    str_to_qual_comb(qual, str, qlen);
#else
    str_to_qual_for_loop(qual, str, qlen);
#endif
}

void qual_to_str(const am_qual_t* qual, char* str, size_t qlen)
{
#if !defined(IN_COMPILER_WE_TRUST)
    qual_to_str_comb(qual, str, qlen);
#else
    qual_to_str_for_loop(qual, str, qlen);
#endif
}

std::string qual_to_str(const std::vector<am_qual_t>& qual)
{
    std::string retq;
    retq.resize(qual.size());
    qual_to_str(qual.data(), retq.data(), qual.size());
    return retq;
}

std::string cigar_arr_to_str(const std::vector<am_cigar_t>& cigar_arr)
{
    return cigar_arr_to_str(cigar_arr.data(), cigar_arr.size());
}

std::string cigar_arr_to_str(const am_cigar_t* cigar_arr, const size_t n)
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
