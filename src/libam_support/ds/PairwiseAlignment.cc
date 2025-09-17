#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include "libam_support/ds/PairwiseAlignment.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"

#include <fmt/format.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
PairwiseAlignment::PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
    std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t pos_on_contig, bool is_plus_strand)
    : aln_query(std::move(aligned_query))
    , aln_ref(std::move(aligned_ref))
    , query(std::move(query))
    , ref(std::move(ref))
    , qual(std::move(qual))
    , read_name(std::move(read_name))
    , contig_name(std::move(contig_name))
    , pos_on_contig(pos_on_contig)
    , is_plus_strand(is_plus_strand)
{
#ifdef CEU_CM_IS_DEBUG
    if (this->aln_ref.length() != this->aln_query.length()) {
        throw PWAException("Length of aligned query and ref inequal!");
    }
    if (this->qual.length() != this->query.length()) {
        throw PWAException("Length of query and qual_ inequal!");
    }
#endif
}

std::vector<am_cigar_t> PairwiseAlignment::generate_cigar_array(const bool use_m) const
{
    const am_cigar_ops_t bam_eq = use_m ? BAM_CMATCH : BAM_CEQUAL;
    const am_cigar_ops_t bam_diff = use_m ? BAM_CMATCH : BAM_CDIFF;

    std::vector<am_cigar_t> cigar;
    am_cigar_ops_t current_cigar_ops = 0;
    am_cigar_ops_t prev_cigar_ops = bam_eq;
    am_cigar_len_t cigar_len = 0;
    const auto ref_len = aln_ref.length();
    for (decltype(aln_ref.length()) i = 0; i < ref_len; i++) {
        if (aln_ref[i] == aln_query[i]) {
            current_cigar_ops = bam_eq;
        } else if (aln_ref[i] == ALN_GAP) {
            current_cigar_ops = BAM_CINS;
        } else if (aln_query[i] == ALN_GAP) {
            current_cigar_ops = BAM_CDEL;
        } else {
            current_cigar_ops = bam_diff;
        }
        if (current_cigar_ops != prev_cigar_ops && cigar_len > 0) {
            cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar_ops));
            cigar_len = 0;
        }
        cigar_len++;
        prev_cigar_ops = current_cigar_ops;
    }
    if (cigar_len != 0) {
        cigar.emplace_back(bam_cigar_gen(cigar_len, prev_cigar_ops));
    }
    return cigar;
}

std::string PairwiseAlignment::serialize() const
{
    return fmt::format(">{}\t{}:{}:{}\n{}\n{}\n{}\n", read_name, contig_name, pos_on_contig, is_plus_strand ? '+' : '-',
        aln_query, aln_ref, qual);
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
    query.erase(std::remove(query.begin(), query.end(), ALN_GAP), query.end());
    ref.erase(std::remove(ref.begin(), ref.end(), ALN_GAP), ref.end());

    return { std::move(read_name), std::move(contig_name), std::move(query), std::move(ref), std::move(qual),
        std::move(aligned_query), std::move(aligned_ref), pos_on_contig, is_plus_strand };
}
[[maybe_unused]] void PairwiseAlignment::serialize(std::ostream& os) const
{
    os << ">" << read_name << "\t" << contig_name << ":" << std::to_string(pos_on_contig) << ":"
       << (is_plus_strand ? '+' : '-') << "\n";
    os << aln_query << "\n";
    os << aln_ref << "\n";
    os << qual << "\n";
}

PWAException::PWAException(const char* msg)
    : msg(msg)
{
}

const char* PWAException::what() const noexcept { return msg; }
} // namespace labw::art_modern
