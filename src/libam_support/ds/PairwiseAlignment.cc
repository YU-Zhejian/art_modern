/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include "libam_support/ds/PairwiseAlignment.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"

#include <fmt/format.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
class PWAException final : public std::runtime_error {
public:
    explicit PWAException(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

PairwiseAlignment::PairwiseAlignment(std::string&& read_name, std::string&& contig_name, std::string&& query,
    std::string&& ref, std::string&& qual_str, std::vector<am_qual_t>&& qual_vec, std::string&& aligned_query,
    std::string&& aligned_ref, hts_pos_t pos_on_contig, bool is_plus_strand)
    : aln_query(std::move(aligned_query))
    , aln_ref(std::move(aligned_ref))
    , query(std::move(query))
    , ref(std::move(ref))
    , qual_str(std::move(qual_str))
    , qual_vec(std::move(qual_vec))
    , read_name(std::move(read_name))
    , contig_name(std::move(contig_name))
    , pos_on_contig(pos_on_contig)
    , is_plus_strand(is_plus_strand)
{
    check_();
}

PairwiseAlignment::PairwiseAlignment(const std::string& read_name, const std::string& contig_name,
    const std::string& query, const std::string& ref, const std::string& qual_str,
    const std::vector<am_qual_t>& qual_vec, const std::string& aligned_query, const std::string& aligned_ref,
    hts_pos_t pos_on_contig, bool is_plus_strand)
    : aln_query(aligned_query)
    , aln_ref(aligned_ref)
    , query(query)
    , ref(ref)
    , qual_str(qual_str)
    , qual_vec(qual_vec)
    , read_name(read_name)
    , contig_name(contig_name)
    , pos_on_contig(pos_on_contig)
    , is_plus_strand(is_plus_strand)
{
    check_();
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

std::string PairwiseAlignment::serialize(const int is_read_1_or_2) const
{
    if (is_read_1_or_2 == 0) {
        return fmt::format(">{}\t{}:{}:{}\n{}\n{}\n{}\n", read_name, contig_name, pos_on_contig,
            (is_plus_strand ? '+' : '-'), aln_query, aln_ref, qual_str);
    }
    return fmt::format(">{}/{}\t{}:{}:{}\n{}\n{}\n{}\n", read_name, is_read_1_or_2, contig_name, pos_on_contig,
        (is_plus_strand ? '+' : '-'), aln_query, aln_ref, qual_str);
}
void PairwiseAlignment::check_() const
{
#ifdef CEU_CM_IS_DEBUG
    if (this->aln_ref.length() != this->aln_query.length()) {
        throw PWAException("Length of aligned query and ref inequal!");
    }
    if (this->qual_str.length() != this->query.length()) {
        throw PWAException("Length of query and qual_ inequal!");
    }
    if (this->qual_vec.size() != this->query.length()) {
        throw PWAException("Length of query and qual_ inequal!");
    }
#endif
}

} // namespace labw::art_modern
