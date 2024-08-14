#include "seq_utils.hh"

#include "art_modern_constants.hh"
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <htslib/sam.h>
#include <map>
#include <sstream>
#include <vector>

using namespace std;
namespace labw {
namespace art_modern {
    const map<char, char>& rev_comp_trans { { 'A', 'T' }, { 'T', 'A' }, { 'C', 'G' }, { 'G', 'C' } };

    std::string qual_to_str(const std::vector<int>& qual)
    {
        string retq;
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
        for (const auto& i : range(0, dna.length(), 1)) {
            try {
                rets[i] = rev_comp_trans.at(dna[i]);
            } catch (std::out_of_range&) {
                // TODO: Make this more user-friendly.
                BOOST_LOG_TRIVIAL(error) << "Invalid character asc(" << (int)dna[i] << ") in dna string";
                throw std::invalid_argument("Invalid character in dna string");
            }
        }
        return rets;
    }
    std::string revcomp(const std::string& dna)
    {
        auto rets = comp(dna);
        reverse(rets.begin(), rets.end());
        return rets;
    }

    std::string cigar_arr_to_str(const vector<uint32_t>& cigar_arr)
    {
        std::ostringstream oss;
        for (auto i = 0; i < cigar_arr.size(); i += 1) {
            oss << (cigar_arr[i] >> BAM_CIGAR_SHIFT);
            oss << BAM_CIGAR_STR[cigar_arr[i] & BAM_CIGAR_MASK];
        }
        return oss.str();
    }

    uint32_t* cigar_arr_to_c(const vector<uint32_t>& cigar_arr)
    {
        auto* c = (uint32_t*)calloc(cigar_arr.size(), sizeof(uint32_t));
        for (auto i = 0; i < cigar_arr.size(); i++) {
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