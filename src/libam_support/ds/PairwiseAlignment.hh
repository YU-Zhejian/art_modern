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

#pragma once

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/hts.h>

#include <array>
#include <exception>
#include <string>
#include <vector>

namespace labw::art_modern {
class PWAException final : public std::exception {
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

    static constexpr int NUM_LINES = 4;

    /**
     *
     * @param read_name Number of read. Will be moved inside.
     * @param contig_name Name of the contig. Will be moved inside.
     * @param query Gapless query sequence. Will be moved inside.
     * @param ref Gapless reference sequence. Will be moved inside.
     * @param qual_str Quality sequence (Added phread offset) whose length should be the same as query. Will be moved
     * inside.
     * @param qual_vec Quality vector whose length should be the same as query. Will be moved inside.
     * @param aligned_query Aligned query sequence with gaps. Will be moved inside.
     * @param aligned_ref Aligned reference sequence with gaps. Will be moved inside.
     * @param pos_on_contig
     * @param is_plus_strand Whether the reference is reverse-complemented.
     */
    PairwiseAlignment(std::string read_name, std::string contig_name, std::string query, std::string ref,
        std::string qual_str, std::vector<am_qual_t> qual_vec, std::string aligned_query, std::string aligned_ref,
        hts_pos_t pos_on_contig, bool is_plus_strand);
    [[maybe_unused]] static PairwiseAlignment deserialize(const std::array<std::string, NUM_LINES>& serialized);
    [[nodiscard]] std::vector<am_cigar_t> generate_cigar_array(bool use_m) const;
    [[nodiscard]] std::string serialize(int is_read_1_or_2 = 0) const;

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
     * Quality sequence whose length should be the same as query.
     */
    const std::string qual_str;
    const std::vector<am_qual_t> qual_vec;
    const std::string read_name;
    const std::string contig_name;
    const hts_pos_t pos_on_contig;
    const bool is_plus_strand;
};

} // namespace labw::art_modern
