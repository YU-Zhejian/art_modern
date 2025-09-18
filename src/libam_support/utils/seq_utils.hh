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

#pragma once
#include "libam_support/Dtypes.hh"

#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

namespace labw::art_modern {
std::string comp(const std::string& dna);
[[maybe_unused]] std::string revcomp(const std::string& dna);
[[maybe_unused]] std::string normalize(const std::string& dna);
void comp_inplace(std::string& dna);
void revcomp_inplace(std::string& dna);
void normalize_inplace(std::string& dna);
std::string qual_to_str(const std::vector<am_qual_t>& qual);
std::string qual_to_str(const am_qual_t* qual, size_t qlen);
std::string qual_to_str_avx2(const am_qual_t* qual, size_t qlen);
std::string qual_to_str_mmx(const am_qual_t* qual, size_t qlen);
std::string qual_to_str_sse2(const am_qual_t* qual, size_t qlen);
std::string qual_to_str_for_loop(const am_qual_t* qual, size_t qlen);
std::string qual_to_str_foreach(const am_qual_t* qual, size_t qlen);
std::string cigar_arr_to_str(const std::vector<am_cigar_t>& cigar_arr);
std::string cigar_arr_to_str(const am_cigar_t* cigar_arr, size_t n);
bool ends_with(const std::string& str, const std::string& suffix);

/**
 * Reverse an arbitrary sequence.
 *
 * @tparam T Sequence type.
 * @param ptr Start of the sequence.
 * @param n Sequence length.
 */
template <typename T> static void reverse(T* ptr, const size_t n)
{
    T tmp;
    for (size_t i = 0; i < n / 2; ++i) {
        tmp = ptr[i];
        ptr[i] = ptr[n - i - 1];
        ptr[n - i - 1] = tmp;
    }
}

template <typename T> std::string vec2str(const std::vector<T>& vec)
{
    std::ostringstream oss;
    std::size_t i = 0;
    oss << "[";
    while (i < vec.size() - 1) {
        oss << std::to_string(vec[i]) + ", ";
        i++;
    }
    oss << std::to_string(vec[i]) << "]";
    return oss.str();
}

} // namespace labw::art_modern
