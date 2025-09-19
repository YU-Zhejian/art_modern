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
#include "art/lib/ArtConstants.hh"
#include "art/lib/Empdist.hh"

#include "libam_support/Constants.hh"

#include <htslib/hts.h>

#include <array>
#include <string>
#include <vector>

namespace labw::art_modern {

struct ArtParams {
    const SIMULATION_MODE art_simulation_mode;
    const ART_LIB_CONST_MODE art_lib_const_mode;
    const bool sep_flag;
    const std::string id;
    const int max_n;
    const int read_len;
    const double pe_frag_dist_mean;
    const double pe_frag_dist_std_dev;
    const std::vector<double> per_base_ins_rate_1;
    const std::vector<double> per_base_del_rate_1;
    const std::vector<double> per_base_ins_rate_2;
    const std::vector<double> per_base_del_rate_2;
    const std::array<double, HIGHEST_QUAL> err_prob;
    const hts_pos_t pe_dist_mean_minus_2_std;
    const Empdist qdist;
};

} // namespace labw::art_modern
