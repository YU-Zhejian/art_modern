#pragma once

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace labw::art_modern {
std::string comp(const std::string& dna);
std::string revcomp(const std::string& dna);
std::string qual_to_str(const std::vector<int>& qual);
std::string qual_to_str(const uint8_t* qual, size_t qlen);
std::string cigar_arr_to_str(const std::vector<uint32_t>& cigar_arr);
std::string normalize(const std::string& dna);

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
    for (int i = 0; i < n / 2; ++i) {
        tmp = ptr[i];
        ptr[i] = ptr[n - i - 1];
        ptr[n - i - 1] = tmp;
    }
}

} // namespace labw::art_modern
