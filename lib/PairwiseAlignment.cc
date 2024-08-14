#include "PairwiseAlignment.hh"
#include "art_modern_constants.hh"
#include "seq_utils.hh"
#include <utility>

namespace labw {
namespace art_modern {

    const char* PWAException::what() const noexcept { return "Alignment is of not equal length!"; }

    PWAException::PWAException(std::string aligned_query, std::string aligned_ref)
        : _aligned_query(std::move(aligned_query))
        , _aligned_ref(std::move(aligned_ref))
    {
    }

    PairwiseAlignment::PairwiseAlignment(std::string read_name, std::string contig_name, std::string query,
        std::string ref, std::string qual, std::string aligned_query, std::string aligned_ref,
        hts_pos_t align_contig_start, bool is_plus_strand)
        : aligned_query(std::move(aligned_query))
        , aligned_ref(std::move(aligned_ref))
        , query(std::move(query))
        , ref(std::move(ref))
        , qual(std::move(qual))
        , read_name(std::move(read_name))
        , contig_name(std::move(contig_name))
        , align_contig_start(align_contig_start)
        , is_plus_strand(is_plus_strand)
    {
        if (this->aligned_ref.length() != this->aligned_query.length()) {
            throw PWAException(this->aligned_query, this->aligned_ref);
        }
    }
    std::vector<uint32_t> PairwiseAlignment::generate_cigar_array(bool use_m) const
    {
        std::vector<uint32_t> cigar;
        uint32_t current_cigar;
        uint32_t prev_cigar;
        uint32_t cigar_len = 0;
        auto ref_len = static_cast<int>(aligned_ref.length());

        auto _range = is_plus_strand ? range(0, ref_len, 1) : range(ref_len - 1, -1, -1);
        for (auto i : _range) {
            if (aligned_ref[i] == aligned_query[i]) {
                current_cigar = use_m ? BAM_CMATCH : BAM_CEQUAL;
            } else if (aligned_ref[i] == ALN_GAP) {
                current_cigar = BAM_CINS;
            } else if (aligned_query[i] == ALN_GAP) {
                current_cigar = BAM_CDEL;
            } else {
                current_cigar = use_m ? BAM_CMATCH : BAM_CDIFF;
            }
            if (current_cigar != prev_cigar && cigar_len > 0) {
                cigar.emplace_back(cigar_len << BAM_CIGAR_SHIFT | prev_cigar);
                cigar_len = 0;
            }
            cigar_len++;
            prev_cigar = current_cigar;
        }
        if (cigar_len != 0) {
            cigar.emplace_back(cigar_len << BAM_CIGAR_SHIFT | prev_cigar);
        }
        return cigar;
    }

} // namespace art_modern
} // namespace labw