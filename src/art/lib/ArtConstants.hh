/**
 * @brief Constants specific to ART.
 *
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
#include "libam_support/Dtypes.hh"

#include <cstdint>

namespace labw::art_modern {

constexpr am_qual_t HIGHEST_QUAL = 80;
constexpr am_qual_count_t MAX_DIST_NUMBER = 1000000;
/**
 * Default insertion rate for read 1.
 */
constexpr double DEFAULT_INS_RATE_1 = 0.00009;
/** Default deletion rate for read 1. */
constexpr double DEFAULT_DEL_RATE_1 = 0.00011;
/**
 * Default insertion rate for read 2.
 */
constexpr double DEFAULT_INS_RATE_2 = 0.00015;
/** Default deletion rate for read 2. */
constexpr double DEFAULT_DEL_RATE_2 = 0.00023;
/** Default maximum number of indels. -1 disables this filter. */
constexpr int DEFAULT_MAX_INDEL = -1;
/** Default maximum number of Ns. 0 allows no N. */
constexpr int DEFAULT_MAX_N = 0;
/** Default batch size. */
constexpr int DEFAULT_BATCH_SIZE = 1U << 14U;

/** Library construction modes. */
enum class ART_LIB_CONST_MODE : std::int8_t {
    /** Single-end mode. */
    SE = 0,
    /** Paired-end mode. */
    PE = 1,
    /** Mate-pair mode. */
    MP = 2
};

constexpr char ART_LIB_CONST_MODE_SE[] = "se";
constexpr char ART_LIB_CONST_MODE_PE[] = "pe";
constexpr char ART_LIB_CONST_MODE_MP[] = "mp";
constexpr char const* ART_LIB_CONST_MODE_STR[]
    = { ART_LIB_CONST_MODE_SE, ART_LIB_CONST_MODE_PE, ART_LIB_CONST_MODE_MP };

/** Program name */
constexpr char ART_PROGRAM_NAME[] = "art_modern";

/**
 * Stop generating reads on a contig if 20% attempts failed.
 */
constexpr double MAX_TRIAL_RATIO_BEFORE_FAIL = 0.2;

} // namespace labw::art_modern