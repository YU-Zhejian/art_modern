#pragma once

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace labw {
namespace art_modern {
    std::string comp(const std::string& dna);
    std::string revcomp(const std::string& dna);
    std::string qual_to_str(const std::vector<int>& qual);
    std::string cigar_arr_to_str(const std::vector<uint32_t>& cigar_arr);
    std::string normalize(const std::string& dna);
    uint32_t* cigar_arr_to_c(const std::vector<uint32_t>& cigar_arr);
} // namespace art_modern
} // namespace labw
