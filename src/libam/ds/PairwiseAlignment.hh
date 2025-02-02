#pragma once

#include "libam/Dtypes.hh"
#include "libam/utils/class_macros_utils.hh"

#include <htslib/hts.h>

#include <array>
#include <exception>
#include <ostream>
#include <string>
#include <vector>

namespace labw::art_modern {
class PWAException : public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override;
    const char* msg;
    explicit PWAException(const char* msg);
};

class PairwiseAlignment {
public:
    DELETE_COPY(PairwiseAlignment)
    DELETE_MOVE(PairwiseAlignment)
    ~PairwiseAlignment() = default;

    static const int NUM_LINES = 4;

    /**
     *
     * @param read_name Number of read. Will be moved inside.
     * @param contig_name Name of the contig. Will be moved inside.
     * @param query Gapless query sequence. Will be moved inside.
     * @param ref Gapless reference sequence. Will be moved inside.
     * @param qual Quality sequence whose length should be the same as query. Will be moved inside.
     * @param aligned_query Aligned query sequence with gaps. Will be moved inside.
     * @param aligned_ref Aligned reference sequence with gaps. Will be moved inside.
     * @param pos_on_contig
     * @param is_plus_strand Whether the reference is reverse-complemented.
     */
    PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
        std::string qual, std::string aligned_query, std::string aligned_ref, hts_pos_t pos_on_contig,
        bool is_plus_strand);
    [[maybe_unused]] static PairwiseAlignment deserialize(const std::array<std::string, NUM_LINES>& serialized);
    [[nodiscard]] std::vector<am_cigar_t> generate_cigar_array(bool use_m) const;
    [[nodiscard]] std::string serialize() const;
    [[maybe_unused]] void serialize(std::ostream& os) const;

    /**
     * Query sequence with gap inserted using -
     */
    const std::string aln_query;
    /**
     * Reference sequence with gap inserted using -
     */
    const std::string aln_ref;
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

} // namespace labw::art_modern
