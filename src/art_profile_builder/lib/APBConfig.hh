/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "art_profile_builder/lib/APBConstants.hh"

#include "libam_support/Dtypes.h"

#include <cstdlib>
#include <string>

namespace labw::art_modern {
/**
 * Art Profile Builder Configuration
 */
class APBConfig {
public:
    const std::string input_file_path;
    const am_read_len_t read_length_1;
    const am_read_len_t read_length_2;
    const std::size_t num_threads;
    const std::size_t num_io_threads;
    const bool is_pe;
    const bool is_ob;
    const std::string output_1_file_path;
    const std::string output_2_file_path;
    const APB_FORMAT format;
    const std::size_t queue_size;
    const am_readnum_t first_n_reads;
};

} // namespace labw::art_modern
