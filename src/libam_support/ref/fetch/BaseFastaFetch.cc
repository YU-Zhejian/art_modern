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

#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include "art_modern_config.h"
#include "libam_support/CExceptionsProxy.hh"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace labw::art_modern {

void BaseFastaFetch::update_sam_header(sam_hdr_t* header) const
{
    for (size_t i = 0; i < seq_lengths_.size(); ++i) {
        auto seq_len_str = std::to_string(seq_lengths_[i]);
        CExceptionsProxy::assert_numeric(
            sam_hdr_add_line(header, "SQ", "SN", seq_names_[i].c_str(), "LN", seq_len_str.c_str(), NULL),
            USED_HTSLIB_NAME, "Failed to add SQ header line for contig '" + seq_names_[i] + "'", false,
            CExceptionsProxy::EXPECTATION::ZERO);
    }
}

hts_pos_t BaseFastaFetch::seq_len(const std::size_t seq_id) const { return seq_lengths_[seq_id]; }

std::string BaseFastaFetch::seq_name(const std::size_t seq_id) const { return seq_names_[seq_id]; }

size_t BaseFastaFetch::num_seqs() const { return seq_lengths_.size(); }

BaseFastaFetch::BaseFastaFetch(std::vector<std::string>&& seq_names, std::vector<hts_pos_t>&& seq_lengths)
    : seq_names_(std::move(seq_names))
    , seq_lengths_(std::move(seq_lengths))
{
}

bool BaseFastaFetch::empty() const { return this->seq_names_.empty(); }
std::string BaseFastaFetch::fetch(const std::size_t seq_id) { return fetch(seq_id, 0, seq_lengths_[seq_id]); }
} // namespace labw::art_modern
