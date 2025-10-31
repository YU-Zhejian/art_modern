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
#include <stdint.h>
#include <stdlib.h>

/**
 * Phred quality score.
 *
 * This score will be involved in +/- operations so it must be signed.
 */
typedef int8_t am_qual_t;

/** Cigar type. Should be 0b00 to 0b11. */
typedef uint8_t am_cigar_type_t;
    /**
 * Cigar operation and length combined.
 */
typedef uint32_t am_cigar_t ;
/**
 * Cigar operation.
 */
typedef uint32_t am_cigar_ops_t;
/**
 * Cigar length.
 */
typedef int32_t am_cigar_len_t;
/** ART quality distribution. */
typedef uint64_t am_qual_count_t;
/**
 * Number of reads.
 *
 * This MUST be a signed integer.
 * Otherwise -1 will be reverted to an extraordinarily large number.
 */
typedef int64_t   am_readnum_t ;
/** Length of files. */
typedef int64_t am_filelen_t;

/** Job ID type. */
typedef size_t am_job_id_t;

