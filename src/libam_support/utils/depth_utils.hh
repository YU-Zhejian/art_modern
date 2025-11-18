/**
 * Copyright 2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "libam_support/Dtypes.h"

#include <cstdlib>
#include <tuple>

namespace labw::art_modern {

std::size_t num_base_positive_pe(
    am_read_len_t read_len_1, am_read_len_t read_len_2, am_readnum_t num_pos_reads, am_readnum_t num_neg_reads);
std::size_t num_base_negative_pe(
    am_read_len_t read_len_1, am_read_len_t read_len_2, am_readnum_t num_pos_reads, am_readnum_t num_neg_reads);
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads_template(
    double cov_pos, double cov_neg, am_readnum_t num_reads_to_reduce);

/**
 * The old algoithm that would work when read lengths are the same.
 */
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads_old(
    std::size_t contig_size, int read_len, double cov_pos, double cov_neg, am_readnum_t num_reads_to_reduce);
} // namespace labw::art_modern
