#pragma once

#include <htslib/hts.h>
#include <string>
#include <vector>

namespace labw::art_modern {
class PWAException : public std::exception {
public:
    const char* what() const noexcept override;
    const char* msg;
    explicit PWAException(const char* msg);
};

class PairwiseAlignment {
public:
    PairwiseAlignment(PairwiseAlignment&& other) = delete;
    PairwiseAlignment(const PairwiseAlignment&) = delete;
    PairwiseAlignment& operator=(PairwiseAlignment&&) = delete;
    PairwiseAlignment& operator=(const PairwiseAlignment&) = delete;
    static const int NUM_LINES = 4;

    /*!
     *
     * @param read_name Number of read.
     * @param contig_name Name of the contig.
     * @param query Gapless query sequence.
     * @param ref Gapless reference sequence.
     * @param qual Quality sequence whose length should be the same as query.
     * @param aligned_query Aligned query sequence with gaps.
     * @param aligned_ref Aligned reference sequence with gaps.
     * @param pos_on_contig
     * @param is_plus_strand Whether the reference is reverse-complemented.
     */
    PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
        std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t pos_on_contig,
        bool is_plus_strand);
    static PairwiseAlignment deserialize(const std::vector<std::string> &serialized);
    std::vector<uint32_t> generate_cigar_array(bool use_m) const;
    std::string serialize() const;

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
    const hts_pos_t pos_on_contig;
    const bool is_plus_strand;

};

} // namespace labw::art_modern // namespace labw