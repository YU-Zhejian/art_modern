#pragma once

#include <string>
#include <vector>

namespace labw {
namespace art_modern {
    std::string comp(const std::string& dna);
    std::string revcomp(const std::string& dna);
    std::string qual_to_str(const std::vector<int>& qual);
} // namespace art_modern
} // namespace labw
