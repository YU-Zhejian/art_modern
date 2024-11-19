#include "seq_utils.hh"

#include "art_modern_constants.hh"
#include "htslib/sam.h"
#include <algorithm>
#include <sstream>
#include <vector>

namespace labw::art_modern {
const char rev_comp_trans_2[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
    24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
    53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
    82, 83, 65, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 };
/*!
 *
 * Generated using:
 *
 * ```python
 * ",".join([str(i) if chr(i) in 'AGCT' else str(ord(chr(i).upper())) if chr(i) in 'agct' else str(ord('N'))for i in
 * range(0, 128)])
 * ```
 */
const char normalization_matrix[] = { 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 65, 78, 67, 78, 78, 78, 71, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 84, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 65, 78, 67, 78, 78, 78, 71, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 84, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78

};

std::string qual_to_str(const uint8_t* qual, const size_t qlen)
{
    std::string retq;
    retq.resize(qlen);
    for (size_t k = 0; k < qlen; k++) {
        retq[k] = (char)(qual[k] + PHRED_OFFSET);
    }
    return retq;
}

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
    for (auto i = 0; i < dna.length(); i++) {
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
    std::for_each(rets.begin(), rets.end(), [](char& c) { c = normalization_matrix[c]; });
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

} // namespace labw::art_modern // namespace labw