#pragma once
#include "libam/Dtypes.hh"

#include <cstddef>
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

/*!
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

} // namespace labw::art_modern
