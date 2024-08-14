#include "seq_utils.hh"

#include <algorithm>
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
            retq[k] = (char)(qual[k] + 33);
        }
        return retq;
    }

    std::string comp(const std::string& dna)
    {
        std::string rets;
        rets.reserve(dna.length());
        for (auto i : dna) {
            rets += rev_comp_trans.at(i);
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
} // namespace art_modern
} // namespace labw