#pragma once

#include <htslib/hts.h>
#include <string>
#include <vector>

namespace labw::art_modern {

class PairwiseAlignment {
public:
    PairwiseAlignment(PairwiseAlignment&& other) = delete;
    PairwiseAlignment(const PairwiseAlignment&) = delete;
    PairwiseAlignment& operator=(PairwiseAlignment&&) = delete;
    PairwiseAlignment& operator=(const PairwiseAlignment&) = delete;

    PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
        std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t align_contig_start,
        bool is_plus_strand);
    std::vector<uint32_t> generate_cigar_array(bool use_m) const;

    /**
     * Query sequence with gap inserted using -
     */
    const std::string aligned_query;
    /**
     * Reference sequence with gap inserted using -
     */
    const std::string aligned_ref;
    /**
     * Query sequence without gap.
     */
    const std::string query;
    /**
     * Reference sequence without gap.
     */
    const std::string ref;
    /**
     * Quality sequence whose elngth should be the same as query.
     */
    const std::string qual;
    const std::string read_name;
    const std::string contig_name;
    const hts_pos_t align_contig_start;
    const bool is_plus_strand;
};

} // namespace labw::art_modern // namespace labw