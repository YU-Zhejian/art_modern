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
#include "art/lib/Rprob.hh"

#include "libam_support/Dtypes.h"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/hts.h>

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace labw::art_modern {

class ArtRead {
public:
    // Disable constructors
    DELETE_MOVE(ArtRead)
    DELETE_COPY(ArtRead)
    DEFAULT_DESTRUCTOR(ArtRead)

    ArtRead(const ArtParams& art_params, std::string contig_name, std::string read_name, bool is_read_1, Rprob& rprob);
    [[nodiscard]] PairwiseAlignment to_pwa();

    int generate_indels();
    // number of deletions <= number of insertions
    int generate_indels_2();

    /**
     * Populate the read while adding insertions and deletions.
     */
    void ref2read(std::string seq_ref, bool is_plus_strand, hts_pos_t pos_on_contig);

    /**
     * Populate aln_read_ and aln_ref_
     */
    void generate_pairwise_aln();

    /**
     * Add point mutations to random bases based on empirical dist of quali scores
     */
    void generate_snv_on_qual();
    [[nodiscard]] bool is_good() const;

private:
    bool is_read_1_;
    am_read_len_t read_len_;
    std::string aln_query_;
    std::string aln_ref_;
    const ArtParams& art_params_;
    std::string contig_name_;
    std::map<int, char, std::less<>> indel_;
    bool is_plus_strand_ = false;
    hts_pos_t pos_on_contig_ = 0;
    std::vector<am_qual_t> qual_;
    std::string read_name_;
    Rprob& rprob_;
    std::string query_;
    std::string ref_;

    [[noreturn]] [[maybe_unused]] void except_() const;
};

} // namespace labw::art_modern
