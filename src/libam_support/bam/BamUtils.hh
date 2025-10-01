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
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/bam/BamTypes.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/sam.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

class BamUtils {
public:
    DELETE_MOVE(BamUtils)
    DELETE_COPY(BamUtils)
    ~BamUtils() = default;
    BamUtils() = delete;

    static void assert_correct_cigar(
        [[maybe_unused]] const PairwiseAlignment& pwa, [[maybe_unused]] const std::vector<am_cigar_t>& cigar);
    static std::string generate_oa_tag(
        const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar, std::int32_t nm_tag);
    static std::pair<std::int32_t, std::string> generate_nm_md_tag(
        const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar);
    static bam1_t_uptr init_uptr();
    static sam_hdr_t* init_header(const BamOptions& sam_options);
    static samFile* open_file(const std::string& filename, const BamOptions& sam_options);
};

} // namespace labw::art_modern