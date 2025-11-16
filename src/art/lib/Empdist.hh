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

#include "art_modern_config.h" // NOLINT

#include "art/lib/Rprob.hh"

#include "libam_support/Dtypes.h"
#include "libam_support/utils/class_macros_utils.hh"

#include "libam_support/ds/GslDiscreteDistribution.hh"

#include <cstddef>
#include <functional>
#include <istream>
#include <map>
#include <string>
#include <vector>

namespace labw::art_modern {

class SlimEmpDistGslDiscrete {
    using dist_map_type = std::map<am_qual_count_t, am_qual_t, std::less<>>;

public:
    DEFAULT_COPY(SlimEmpDistGslDiscrete)
    DEFAULT_MOVE(SlimEmpDistGslDiscrete)
    DEFAULT_DESTRUCTOR(SlimEmpDistGslDiscrete)
    explicit SlimEmpDistGslDiscrete(const dist_map_type& dist)
    {
        std::vector<double> count;
        for (const auto& [this_count, this_qual] : dist) {
            qual_.emplace_back(this_qual);
            count.emplace_back(static_cast<double>(this_count));
        }
        std::vector<double> init_list;
        double prev = 0;
        for (const auto i : count) {
            init_list.emplace_back(i - prev);
            prev = i;
        }
        rd_ = GslDiscreteDistribution(init_list);
    }

    [[nodiscard]] am_qual_t gen_qual(const double u) const { return qual_[rd_(u)]; }

private:
    std::vector<am_qual_t> qual_;
    GslDiscreteDistribution<double> rd_;
};

class Empdist {
public:
    using dist_map_type = std::map<am_qual_count_t, am_qual_t, std::less<>>;
    using dist_type = std::vector<dist_map_type>;

private:
    dist_type qual_dist_first;
    dist_type qual_dist_second;

    dist_type a_qual_dist_first;
    dist_type t_qual_dist_first;
    dist_type g_qual_dist_first;
    dist_type c_qual_dist_first;

    dist_type a_qual_dist_second;
    dist_type t_qual_dist_second;
    dist_type g_qual_dist_second;
    dist_type c_qual_dist_second;
    using dist_idx_type = std::vector<SlimEmpDistGslDiscrete>;
    dist_idx_type qual_dist_first_idx;
    dist_idx_type qual_dist_second_idx;
    dist_idx_type a_qual_dist_first_idx;
    dist_idx_type t_qual_dist_first_idx;
    dist_idx_type g_qual_dist_first_idx;
    dist_idx_type c_qual_dist_first_idx;
    dist_idx_type a_qual_dist_second_idx;
    dist_idx_type t_qual_dist_second_idx;
    dist_idx_type g_qual_dist_second_idx;
    dist_idx_type c_qual_dist_second_idx;

public:
    Empdist(const std::string& emp_filename_1, const std::string& emp_filename_2, bool sep_qual, bool is_pe);
    Empdist(const std::string& builtin_profile_name, bool sep_qual, bool is_pe, bool silence = false);

    void set_read_length(am_read_len_t read_len_1, am_read_len_t read_len_2);
    void shift_all_emp(am_qual_t q_shift_1, am_qual_t q_shift_2, am_qual_t min_qual, am_qual_t max_qual);
    /** Prepare the WALKER indices **/
    void index();
    /** Log the loaded profile **/
    void log() const;

    void get_read_qual(std::vector<am_qual_t>& qual, Rprob& rprob, bool first = true) const;
    void get_read_qual_sep_1(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const;
    void get_read_qual_sep_2(std::vector<am_qual_t>& qual, const std::string& seq, Rprob& rprob) const;
    [[nodiscard]] am_read_len_t get_read_1_max_length() const;
    [[nodiscard]] am_read_len_t get_read_2_max_length() const;

private:
    void read_emp_dist_(const std::string& infile, bool is_first);
    void read_emp_dist_(std::istream& input, bool is_first);
    void validate_() const;
    const bool sep_qual_;
    const bool is_pe_;
    am_read_len_t read_len_1_ { 0 };
    am_read_len_t read_len_2_ { 0 };
    const bool silence_;
};

} // namespace labw::art_modern
