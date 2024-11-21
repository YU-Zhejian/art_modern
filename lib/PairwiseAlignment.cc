#include "PairwiseAlignment.hh"
#include "art_modern_config.h" // For CEU_CM_IS_DEBUG
#include "art_modern_constants.hh"
#include <htslib/sam.h>
#include <utility>

namespace labw::art_modern {
PairwiseAlignment::PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
    std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t align_contig_start,
    bool is_plus_strand)
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
#ifdef CEU_CM_IS_DEBUG
    if (this->aligned_ref.length() != this->aligned_query.length()) {
        throw PWAException("Length of aligned query and ref inequal!");
    }
    if (this->qual.length() != this->query.length()) {
        throw PWAException("Length of query and qual inequal!");
    }
#endif
}

std::vector<uint32_t> PairwiseAlignment::generate_cigar_array(const bool use_m) const
{
    std::vector<uint32_t> cigar;
    uint32_t current_cigar = BAM_CMATCH;
    uint32_t prev_cigar = BAM_CMATCH;
    uint32_t cigar_len = 0;
    const auto ref_len = aligned_ref.length();
    if (is_plus_strand) {
        for (auto i = 0; i < ref_len; i++) {
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
                cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar));
                cigar_len = 0;
            }
            cigar_len++;
            prev_cigar = current_cigar;
        }
    } else {
        for (auto i = static_cast<int>(ref_len) - 1; i > -1; i--) {
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
                cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar));
                cigar_len = 0;
            }
            cigar_len++;
            prev_cigar = current_cigar;
        }
    }
    if (cigar_len != 0) {
        cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar));
    }
    return cigar;
}

PWAException::PWAException(const char* msg)
    : msg(msg)
{
}

const char* PWAException::what() const noexcept { return msg; }
} // namespace labw::art_modern // namespace labw