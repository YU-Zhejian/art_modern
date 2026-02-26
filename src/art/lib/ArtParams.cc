/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
 * Copyright 2024-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
#include "art/lib/ArtParams.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/Empdist.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/arithmetic_utils.hh"

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

namespace {

    std::array<double, HIGHEST_QUAL> gen_err_prob_()
    {
        std::array<double, HIGHEST_QUAL> tmp_err_prob {};
        for (int i = 0; i < HIGHEST_QUAL; i++) {
            tmp_err_prob[i] = std::pow(10, -i / 10.0);
        }
        return tmp_err_prob;
    }
} // namespace
ArtParams::ArtParams(const SIMULATION_MODE art_simulation_mode, const ART_LIB_CONST_MODE art_lib_const_mode,
    const bool sep_flag, std::string&& id, const int max_n, const am_read_len_t read_len_1,
    const am_read_len_t read_len_2, const double pe_frag_dist_mean, const double pe_frag_dist_std_dev,
    std::vector<double>&& per_base_ins_rate_1, std::vector<double>&& per_base_del_rate_1,
    std::vector<double>&& per_base_ins_rate_2, std::vector<double>&& per_base_del_rate_2, Empdist&& qdist,
    const std::size_t job_pool_reporting_interval_seconds,
    const std::size_t art_job_executor_reporting_interval_seconds, const am_rand_seed_t seed)
    : art_simulation_mode(art_simulation_mode)
    , art_lib_const_mode(art_lib_const_mode)
    , sep_flag(sep_flag)
    , id(std::move(id))
    , max_n(max_n)
    , read_len_1(read_len_1)
    , read_len_2(read_len_2)
    , pe_frag_dist_mean(pe_frag_dist_mean)
    , pe_frag_dist_std_dev(pe_frag_dist_std_dev)
    , per_base_ins_rate_1(std::move(per_base_ins_rate_1))
    , per_base_del_rate_1(std::move(per_base_del_rate_1))
    , per_base_ins_rate_2(std::move(per_base_ins_rate_2))
    , per_base_del_rate_2(std::move(per_base_del_rate_2))
    , qdist(std::move(qdist))
    , job_pool_reporting_interval_seconds(job_pool_reporting_interval_seconds)
    , art_job_executor_reporting_interval_seconds(art_job_executor_reporting_interval_seconds)
    , err_prob(gen_err_prob_())
    , pe_dist_mean_minus_2_std(static_cast<am_read_len_t>(pe_frag_dist_mean - 2 * pe_frag_dist_std_dev))
    , seed(seed)
{
    if (art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
        contig_len_threshold = read_len_1;
        return;
    }
    if (art_lib_const_mode == ART_LIB_CONST_MODE::PE || art_lib_const_mode == ART_LIB_CONST_MODE::MP) {
        if (art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
            contig_len_threshold = am_max(read_len_1, read_len_2);
        } else {
            contig_len_threshold = am_max(pe_dist_mean_minus_2_std, am_max(read_len_1, read_len_2));
        }
        return;
    }
}

} // namespace labw::art_modern
