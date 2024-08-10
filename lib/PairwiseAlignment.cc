//
// Created by yuzj on 3/22/24.
//

#include "PairwiseAlignment.hh"
#include "art_modern_constants.hh"
#include "misc.hh"
#include <boost/format.hpp>
#include <sstream>
#include <utility>

namespace labw {
namespace art_modern {

    const char* PWAException::what() const noexcept
    {
        return (boost::format("Alignment %s -> %s is of not equal length!") % _aligned_query % _aligned_ref).str().c_str();
    }

    PWAException::PWAException(std::string aligned_query, std::string aligned_ref)
        : _aligned_query(std::move(aligned_query))
        , _aligned_ref(std::move(aligned_ref))
    {
    }

    std::string PairwiseAlignment::generate_cigar(bool is_reverse, bool use_m) const
    {
        std::ostringstream cigar;

        char current_cigar;
        char prev_cigar;
        int cigar_len = 0;
        auto ref_len = static_cast<int>(_aligned_ref.length());

        auto _range = is_reverse ? labw::art_modern::range(ref_len - 1, -1, -1) : labw::art_modern::range(0, ref_len, 1);
        for (auto i : _range) {
            if (_aligned_ref[i] == _aligned_query[i]) {
                current_cigar = use_m ? ALN_MATCH : SEQ_MATCH;
            } else if (_aligned_ref[i] == ALN_GAP) {
                current_cigar = INSERTION;
            } else if (_aligned_query[i] == ALN_GAP) {
                current_cigar = DELETION;
            } else {
                current_cigar = use_m ? ALN_MATCH : SEQ_MISMATCH;
            }
            if (current_cigar != prev_cigar && cigar_len > 0) {
                cigar << cigar_len << prev_cigar;
                cigar_len = 0;
            }
            cigar_len++;
            prev_cigar = current_cigar;
        }
        if (cigar_len != 0) {
            cigar << cigar_len << prev_cigar;
        }
        return cigar.str();
    }

    PairwiseAlignment::PairwiseAlignment(const std::string& aligned_query, const std::string& aligned_ref)
        : _aligned_query(aligned_query)
        , _aligned_ref(aligned_ref)
    {
        if (_aligned_ref.length() != _aligned_query.length()) {
            throw PWAException(aligned_query, aligned_ref);
        }
    }
}
}