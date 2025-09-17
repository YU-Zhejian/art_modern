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

#include "art/lib/ArtContig.hh"

#include "art/lib/ArtRead.hh"
#include "art/lib/Rprob.hh"

#include "libam_support/Constants.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <htslib/hts.h>

#include <cstddef>
#include <memory>
#include <utility>

namespace labw::art_modern {

/**
 * SE:@code
                |----------->
             ------------------------------------
 OR
             ------------------------------------
                               <-----------|
 * @endcode
 * @param is_plus_strand
 * @param read_1
 */
void ArtContig::generate_read_se(const bool is_plus_strand, ArtRead& read_1)
{
    const auto pos_1 = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        ? 0
        : rprob_.randint(0, static_cast<int>(valid_region_) + 1);
    auto slen_1 = read_1.generate_indels(true);
    // ensure get a fixed read length
    if (pos_1 + art_params_.read_len - slen_1 > seq_size) {
        slen_1 = read_1.generate_indels_2(true);
    }
    auto seq_ref = fasta_fetch_->fetch(seq_id_, pos_1, pos_1 + art_params_.read_len - slen_1);
    read_1.ref2read(std::move(seq_ref), is_plus_strand, pos_1);
}

/**
 * PE: @code
                |----------->
             ------------------------------------
                               <-----------|

 * @endcode
 * MP: @code
                <-----------|
             ------------------------------------
                               |----------->

 * @endcode
 *
 * @param is_plus_strand
 * @param is_mp
 * @param read_1
 * @param read_2
 */
void ArtContig::generate_read_pe(const bool is_plus_strand, const bool is_mp, ArtRead& read_1, ArtRead& read_2)
{
    const hts_pos_t fragment_len = generate_fragment_length();
    const hts_pos_t fragment_start = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        ? 0
        : rprob_.randint(0, static_cast<int>(seq_size - fragment_len) + 1);
    const hts_pos_t fragment_end = fragment_start + fragment_len;

    const hts_pos_t pos_1 = is_mp == is_plus_strand ? fragment_end - art_params_.read_len : fragment_start;
    const hts_pos_t pos_2 = is_mp == is_plus_strand ? fragment_start : fragment_end - art_params_.read_len;

    int slen_1 = read_1.generate_indels(true);
    int slen_2 = read_2.generate_indels(false);

    // ensure get a fixed read length
    if (pos_1 + art_params_.read_len - slen_1 > seq_size) {
        slen_1 = read_1.generate_indels_2(true);
    }
    if (pos_2 + art_params_.read_len - slen_2 > seq_size) {
        slen_2 = read_2.generate_indels_2(false);
    }
    auto seq_ref_1 = fasta_fetch_->fetch(seq_id_, pos_1, pos_1 + art_params_.read_len - slen_1);
    auto seq_ref_2 = fasta_fetch_->fetch(seq_id_, pos_2, pos_2 + art_params_.read_len - slen_2);

    read_1.ref2read(std::move(seq_ref_1), is_plus_strand, pos_1);
    read_2.ref2read(std::move(seq_ref_2), !is_plus_strand, pos_2);
}

ArtContig::ArtContig(
    const std::shared_ptr<BaseFastaFetch>& fasta_fetch, const size_t seq_id, const ArtParams& art_params, Rprob& rprob)
    : seq_name(fasta_fetch->seq_name(seq_id))
    , seq_size(fasta_fetch->seq_len(seq_id))
    , art_params_(art_params)
    , fasta_fetch_(fasta_fetch)
    , rprob_(rprob)
    , seq_id_(seq_id)
    , valid_region_(fasta_fetch->seq_len(seq_id) - art_params.read_len)
{
}
hts_pos_t ArtContig::generate_fragment_length() const
{
    hts_pos_t fragment_len = 0;
    if (art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
        || art_params_.pe_dist_mean_minus_2_std > seq_size) {
        // when reference length < pe_frag_dist_mean-2*std, fragment_len sets to
        // be reference length
        fragment_len = seq_size;
    } else {
        fragment_len = 0;
        while (fragment_len < art_params_.read_len || fragment_len > seq_size) {
            fragment_len = rprob_.insertion_length();
        }
    }
    return fragment_len;
}

} // namespace labw::art_modern
