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

#include "art_modern_config.h"
#include "libam_support/Constants.hh"

#include <string>

namespace labw::art_modern {
class BamOptions {
public:
    static constexpr char ALLOWED_COMPRESSION_LEVELS[] { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'u' };
    /**
     * Format version. Accepted format: `/^[0-9]+\.[0-9]+$`.
     */
    std::string HD_VN = "1.4";
    /**
     * Sorting order of alignments. Valid values: `unknown` (default), `unsorted`, `queryname` and `coordinate`.
     */
    std::string HD_SO = "unsorted";

    /**
     * Program record identifier.
     */
    std::string PG_ID = "01";
    /**
     * Program name.
     */
    std::string PG_PN = "art_modern";
    /**
     * Command line.
     */
    std::string PG_CL;
    /**
     * Program version.
     */
    std::string PG_VN = std::string("ART-") + ART_VERSION + "-ART_MODERN-" + ART_MODERN_VERSION;

    /**
     * Use `M` instead of `X` or `=` for matching.
     */
    bool use_m = false;

    /**
     * If `false`, will write SAM instead.
     */
    bool write_bam = true;
    /**
     * Number of threads used by htslib.
     */
    int hts_io_threads = 1;

    /**
     * Compression level `[u0-9]`.
     */
    char compress_level = '4';
};
} // namespace labw::art_modern
