#include "seq_utils.hh"

#include "art_modern_constants.hh"
#include <algorithm>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/log/trivial.hpp>
#include <htslib/sam.h>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace labw {
namespace art_modern {
    const char rev_comp_trans_2[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
        23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
        51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73, 74, 75, 76, 77, 78,
        79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104,
        105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
        127 };

    std::string qual_to_str(const std::vector<int>& qual)
    {
        std::string retq;
        retq.resize(qual.size());
        for (size_t k = 0; k < qual.size(); k++) {
            retq[k] = (char)(qual[k] + PHRED_OFFSET);
        }
        return retq;
    }

    std::string comp(const std::string& dna)
    {
        std::string rets;
        rets.resize(dna.length());
        for (const auto& i : range(0, static_cast<int>(dna.length()), 1)) {
            rets[i] = rev_comp_trans_2[dna[i]];
        }
        return rets;
    }
    std::string revcomp(const std::string& dna)
    {
        auto rets = comp(dna);
        reverse(rets.begin(), rets.end());
        return rets;
    }

    std::string normalize(const std::string& dna)
    {
        std::string rets = dna;
        boost::algorithm::to_upper(rets);
        for (const auto& i : range(0, static_cast<int>(rets.length()), 1)) {
            if (rets[i] != 'A' && rets[i] != 'T' && rets[i] != 'C' && rets[i] != 'G') {
                rets[i] = 'N';
            }
        }
        return rets;
    }

    std::string cigar_arr_to_str(const std::vector<uint32_t>& cigar_arr)
    {
        std::ostringstream oss;
        for (size_t i = 0; i < cigar_arr.size(); i += 1) {
            oss << (cigar_arr[i] >> BAM_CIGAR_SHIFT);
            oss << BAM_CIGAR_STR[cigar_arr[i] & BAM_CIGAR_MASK];
        }
        return oss.str();
    }

    uint32_t* cigar_arr_to_c(const std::vector<uint32_t>& cigar_arr)
    {
        auto* c = (uint32_t*)calloc(cigar_arr.size(), sizeof(uint32_t));
        for (size_t i = 0; i < cigar_arr.size(); i++) {
            c[i] = cigar_arr[i];
        }
        return c;
    }

    std::vector<int> range(int start, int stop, int step)
    {
        std::vector<int> result;
        for (int i = start; step > 0 ? i < stop : i > stop; i += step) {
            result.push_back(i);
        }
        return result;
    }
} // namespace art_modern
} // namespace labw