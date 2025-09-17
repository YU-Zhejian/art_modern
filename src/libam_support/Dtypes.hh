/**
 *@brief Basic data types.
 *
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
#include <cstdint>

namespace labw::art_modern {
/**
 * Phred quality score.
 *
 * This score will be involved in +/- operations so it must be signed.
 */
using am_qual_t = int8_t;
/** Cigar type. Should be 0b00 to 0b11. */
using am_cigar_type_t = uint8_t;
/**
 * Cigar operation and length combined.
 */
using am_cigar_t = uint32_t;
/**
 * Cigar operation.
 */
using am_cigar_ops_t = uint32_t;
/**
 * Cigar length.
 */
using am_cigar_len_t = int32_t;
/** ART quality distribution. */
using am_qual_count_t = uint64_t;
/**
 * Number of reads.
 *
 * This MUST be a signed integer.
 * Otherwise -1 will be reverted to an extraordinarily large number.
 */
using am_readnum_t = int64_t;
/** Length of files. */
using am_filelen_t = int64_t;
} // namespace labw::art_modern
