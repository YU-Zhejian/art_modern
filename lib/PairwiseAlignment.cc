#include "PairwiseAlignment.hh"
#include "art_modern_config.h" // For CEU_CM_IS_DEBUG
#include "art_modern_constants.hh"
#include "utils/seq_utils.hh"

#include <boost/algorithm/string/erase.hpp>
#include <htslib/sam.h>
#include <sstream>
#include <utility>

namespace labw::art_modern {
PairwiseAlignment::PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
    std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t pos_on_contig, bool is_plus_strand)
    : aligned_query(std::move(aligned_query))
    , aligned_ref(std::move(aligned_ref))
    , query(std::move(query))
    , ref(std::move(ref))
    , qual(std::move(qual))
    , read_name(std::move(read_name))
    , contig_name(std::move(contig_name))
    , pos_on_contig(pos_on_contig)
    , is_plus_strand(is_plus_strand)
{
#ifdef CEU_CM_IS_DEBUG
    if (this->aligned_ref.length() != this->aligned_query.length()) {
        throw PWAException("Length of aligned query and ref inequal!");
    }
    if (this->qual.length() != this->query.length()) {
        throw PWAException("Length of query and qual_ inequal!");
    }
#endif
}

std::vector<uint32_t> PairwiseAlignment::generate_cigar_array(const bool use_m) const
{
    const uint32_t bam_eq = use_m ? BAM_CMATCH : BAM_CEQUAL;
    const uint32_t bam_diff= use_m ? BAM_CMATCH : BAM_CDIFF;

    std::vector<uint32_t> cigar;
    uint32_t current_cigar;
    uint32_t prev_cigar = bam_eq;
    uint32_t cigar_len = 0;
    const auto ref_len = aligned_ref.length();
    for (decltype(aligned_ref.length()) i = 0; i < ref_len; i++) {
        if (aligned_ref[i] == aligned_query[i]) {
            current_cigar = bam_eq;
        } else if (aligned_ref[i] == ALN_GAP) {
            current_cigar = BAM_CINS;
        } else if (aligned_query[i] == ALN_GAP) {
            current_cigar = BAM_CDEL;
        } else {
            current_cigar = bam_diff;
        }
        if (current_cigar != prev_cigar && cigar_len > 0) {
            cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar));
            cigar_len = 0;
        }
        cigar_len++;
        prev_cigar = current_cigar;
    }
    if (cigar_len != 0) {
        cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar));
    }
    return cigar;
}

std::string PairwiseAlignment::serialize() const
{
    std::ostringstream os;
    os << ">" << read_name << "\t" << contig_name << ":" << std::to_string(pos_on_contig) << ":"
       << (is_plus_strand ? '+' : '-') << "\n";
    os << aligned_query << "\n";
    os << aligned_ref << "\n";
    os << qual << "\n";
    return os.str();
}

PairwiseAlignment PairwiseAlignment::deserialize(const std::vector<std::string>& serialized)
{
    const auto sep_pos = serialized[0].find('\t');
    const std::string& read_name_ = serialized[0].substr(1, sep_pos - 1);
    const std::string& coordinate_ = serialized[0].substr(sep_pos + 1);
    std::string token;
    std::istringstream iss(coordinate_);
    std::getline(iss, token, ':');
    const std::string& contig_name_ = token;
    std::getline(iss, token, ':');
    const hts_pos_t pos_on_contig_ = std::stol(token);
    std::getline(iss, token, ':');
    const bool is_plus_strand_ = token[0] == '+';

    const std::string& aligned_query_ = serialized[1];
    const std::string& aligned_ref_ = serialized[2];
    std::string query_ = serialized[1];
    std::string ref_ = serialized[2];
    const std::string& qual_ = serialized[3];
    boost::algorithm::erase_all(query_, ALN_GAP_STR);
    boost::algorithm::erase_all(ref_, ALN_GAP_STR);
    return { read_name_, contig_name_, query_, ref_, qual_, aligned_query_, aligned_ref_, pos_on_contig_,
        is_plus_strand_ };
}

PWAException::PWAException(const char* msg)
    : msg(msg)
{
}

const char* PWAException::what() const noexcept { return msg; }
} // namespace labw::art_modern // namespace labw