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

#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace labw::art_modern {

enum class FastaFetchType : std::uint8_t {
    IN_MEMORY_FASTA_FETCH,
    FAIDX_FETCH,
};

/**
 * Random access FASTA file.
 */
class BaseFastaFetch {
public:
    BaseFastaFetch(const BaseFastaFetch& other)
        : BaseFastaFetch(other.seq_names_, other.seq_lengths_)
    {
    }
    BaseFastaFetch& operator=(const BaseFastaFetch&) = delete;

    DELETE_MOVE(BaseFastaFetch)

    BaseFastaFetch() = default;

    /**
     * Default destructor.
     */
    virtual ~BaseFastaFetch() = default;

    explicit BaseFastaFetch(const std::tuple<std::vector<std::string>, std::vector<hts_pos_t>>& seq_names_lengths);
    BaseFastaFetch(std::vector<std::string>&& seq_names, std::vector<hts_pos_t>&& seq_lengths);
    BaseFastaFetch(const std::vector<std::string>& seq_names, const std::vector<hts_pos_t>& seq_lengths);

    /**
     * This method is thread-safe since mutex is used for non-thread-safe implementations.
     *
     * @param seq_id Contig name.
     * @param start 0-based inclusive start point.
     * @param end 0-based exclusive end point.
     * @return Fetched sequence.
     */
    virtual std::string fetch(std::size_t seq_id, hts_pos_t start, hts_pos_t end) = 0;

    /**
     * This method is thread-safe since mutex is used for non-thread-safe implementations.
     * Fetch the entire contig.
     *
     * @param seq_id Contig name.
     * @return Fetched sequence.
     */
    virtual std::string fetch(std::size_t seq_id);

    void update_sam_header(sam_hdr_t* header) const;

    /**
     * Get the length of the desired sequence.
     * This operation should take $O(1)$ complexity.
     *
     * @param seq_id As described.
     * @return As described.
     */
    [[nodiscard]] hts_pos_t seq_len(std::size_t seq_id) const;

    /**
     * Get the name of the desired sequence.
     * @param seq_id As described.
     * @return As described.
     */
    [[nodiscard]] std::string seq_name(std::size_t seq_id) const;

    /**
     * Get number of sequences inside.
     *
     * @return As described.
     */
    [[nodiscard]] size_t num_seqs() const;

    [[nodiscard]] bool empty() const;

protected:
    const std::vector<std::string> seq_names_;
    const std::vector<hts_pos_t> seq_lengths_;
};
} // namespace labw::art_modern
