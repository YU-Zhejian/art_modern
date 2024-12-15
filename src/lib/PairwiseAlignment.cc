#include "PairwiseAlignment.hh"
#include "art_modern_config.h" // For CEU_CM_IS_DEBUG
#include "art_modern_constants.hh"

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

std::vector<am_cigar_t> PairwiseAlignment::generate_cigar_array(const bool use_m) const
{
    const am_cigar_t bam_eq = use_m ? BAM_CMATCH : BAM_CEQUAL;
    const am_cigar_t bam_diff = use_m ? BAM_CMATCH : BAM_CDIFF;

    std::vector<am_cigar_t> cigar;
    am_cigar_t current_cigar;
    am_cigar_t prev_cigar = bam_eq;
    am_cigar_t cigar_len = 0;
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
    const auto pos_on_contig_s = std::to_string(pos_on_contig);
    const std::size_t buff_len = read_name.length() + contig_name.length() + pos_on_contig_s.length()
        + aligned_query.length() + aligned_ref.length() + qual.length() + 9;
    std::string buff;
    buff.resize(buff_len);
    buff[0] = 0;
    std::snprintf(buff.data(), buff_len + 1, ">%s\t%s:%s:%c\n%s\n%s\n%s\n", read_name.c_str(), contig_name.c_str(),
        pos_on_contig_s.c_str(), is_plus_strand ? '+' : '-', aligned_query.c_str(), aligned_ref.c_str(), qual.c_str());

    return buff;
}

[[maybe_unused]] PairwiseAlignment PairwiseAlignment::deserialize(const std::array<std::string, NUM_LINES>& serialized)
{
    const auto sep_pos = serialized[0].find('\t');
    std::string read_name = serialized[0].substr(1, sep_pos - 1);
    const std::string& coordinate = serialized[0].substr(sep_pos + 1);
    std::string token;
    std::istringstream iss(coordinate);
    std::getline(iss, token, ':');
    std::string contig_name = std::move(token);
    std::getline(iss, token, ':');
    const hts_pos_t pos_on_contig = std::stol(token);
    std::getline(iss, token, ':');
    const bool is_plus_strand = token[0] == '+';

    std::string aligned_query = serialized[1];
    std::string aligned_ref = serialized[2];
    std::string query = serialized[1];
    std::string ref = serialized[2];
    std::string qual = serialized[3];
    boost::algorithm::erase_all(query, ALN_GAP_STR);
    boost::algorithm::erase_all(ref, ALN_GAP_STR);
    return { std::move(read_name), std::move(contig_name), std::move(query), std::move(ref), std::move(qual),
        std::move(aligned_query), std::move(aligned_ref), pos_on_contig, is_plus_strand };
}
[[maybe_unused]] void PairwiseAlignment::serialize(std::ostream& os) const
{
    os << ">" << read_name << "\t" << contig_name << ":" << std::to_string(pos_on_contig) << ":"
       << (is_plus_strand ? '+' : '-') << "\n";
    os << aligned_query << "\n";
    os << aligned_ref << "\n";
    os << qual << "\n";
}

PWAException::PWAException(const char* msg)
    : msg(msg)
{
}

const char* PWAException::what() const noexcept { return msg; }
} // namespace labw::art_modern // namespace labw