/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
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

#include "art/lib/ArtParams.hh"
#include "art/lib/ArtRead.hh"
#include "art/lib/Rprob.hh"

#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <htslib/hts.h>

#include <cstdlib>
#include <memory>
#include <string>

namespace labw::art_modern {

class ArtContig {

public:
    ArtContig(
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch, size_t seq_id, const ArtParams& art_params, Rprob& rprob);

    void generate_read_se(bool is_plus_strand, ArtRead& read_1);
    void generate_read_pe(bool is_plus_strand, bool is_mp, ArtRead& read_1, ArtRead& read_2);

    /** Name of the contig */
    const std::string seq_name;
    /**Size of the contig */
    const hts_pos_t seq_size;

private:
    [[nodiscard]] hts_pos_t generate_fragment_length() const;
    const ArtParams& art_params_;
    const std::shared_ptr<BaseFastaFetch>& fasta_fetch_;
    Rprob& rprob_;
    const std::size_t seq_id_;
    const hts_pos_t valid_region_;
};

} // namespace labw::art_modern