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

#include "art_modern_config.h" // NOLINT: For USE_STL_RANDOM, etc.

#include "libam_support/Dtypes.h"
#include "libam_support/utils/class_macros_utils.hh"

#if defined(USE_STL_RANDOM)
// Do nothing, included below
#elif defined(USE_BOOST_RANDOM)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#elif defined(USE_ONEMKL_RANDOM)
#include <mkl.h> // NOLINT
#elif defined(USE_PCG_RANDOM)
#include "libam_support/ds/pcg_32_c.hh"
#define PCG_CLASS_NAME pcg32_c
#elif defined(USE_SYSTEM_PCG_RANDOM)
#include <pcg_random.hpp>
#define PCG_CLASS_NAME pcg32_fast
#else
#error                                                                                                                 \
    "Define USE_STL_RANDOM, USE_BOOST_RANDOM, USE_ONEMKL_RANDOM, USE_PCG_RANDOM, USE_SYSTEM_PCG_RANDOM for random generators!"
#endif

#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_SYSTEM_PCG_RANDOM)
#define USE_STL_LIKE_RANDOM
#endif

#include <htslib/hts.h>

// Include C++ stdlibs
#if defined(USE_ONEMKL_RANDOM)
#include <array> // Cache uses std::array
#endif
#include <cstddef>
#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM) || defined(USE_SYSTEM_PCG_RANDOM)
#include <random>
#endif
#include <vector>

#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM) || defined(USE_SYSTEM_PCG_RANDOM)
#define REAL_DIST std::uniform_real_distribution
#define NORM_DIST std::normal_distribution
#define INT_DIST std::uniform_int_distribution
#elif defined(USE_BOOST_RANDOM)
#define REAL_DIST boost::random::uniform_01
#define NORM_DIST boost::random::normal_distribution
#define INT_DIST boost::random::uniform_int_distribution
#endif

namespace labw::art_modern {

class Rprob {
public:
    Rprob(double pe_frag_dist_mean, double pe_frag_dist_std_dev, am_read_len_t read_len_1, am_read_len_t read_len_2);
    DELETE_COPY(Rprob)
    DELETE_MOVE(Rprob)
    ~Rprob();

    /**
     * Populate tmp_probs_ with n random values using r_prob().
     *
     * n must be less than or equal to read_length_.
     */
    void r_probs(std::size_t n);
    /**
     * Generate an insertion length based on Gaussian distribution.
     * @return
     */
    hts_pos_t insertion_length();
    /**
     * Generate one of A, C, G, T.
     */
    char rand_base();
    am_qual_t rand_quality_less_than_10();
    int rand_pos_on_read(bool is_read1);
    int rand_pos_on_read_not_head_and_tail(bool is_read1);
    /** Slow random integer function **/
    int randint(int min, int max);
    std::vector<double> tmp_probs_;

private:
    void public_init_();

#if defined(USE_ONEMKL_RANDOM)
    // OneMKL bulk random number generation cache
    constexpr static std::size_t CACHE_SIZE_ = 4096;

    std::array<int, CACHE_SIZE_> cached_rand_pos_on_read_1_ {};
    std::size_t cached_rand_pos_on_read_1_index_ = 0;

    std::array<int, CACHE_SIZE_> cached_rand_pos_on_read_2_ {};
    std::size_t cached_rand_pos_on_read_2_index_ = 0;

    std::array<int, CACHE_SIZE_> cached_rand_pos_on_read_1_not_head_and_tail_ {};
    std::size_t cached_rand_pos_on_read_1_not_head_and_tail_index_ = 0;

    std::array<int, CACHE_SIZE_> cached_rand_pos_on_read_2_not_head_and_tail_ {};
    std::size_t cached_rand_pos_on_read_2_not_head_and_tail_index_ = 0;

    std::array<double, CACHE_SIZE_> cached_insertion_lengths_ {};
    std::size_t cached_insertion_lengths_index_ = 0;

    std::array<int, CACHE_SIZE_> cached_rand_base_indices_ {};
    std::size_t cached_rand_base_indices_index_ = 0;

    std::array<am_qual_t, CACHE_SIZE_> cached_rand_quality_less_than_10_ {};
    std::size_t cached_rand_quality_less_than_10_index_ = 0;
#endif

#if defined(USE_STL_LIKE_RANDOM)
#if defined(USE_STL_RANDOM)
    std::mt19937 gen_;
#elif defined(USE_PCG_RANDOM) || defined(USE_SYSTEM_PCG_RANDOM)
    PCG_CLASS_NAME gen_;
#else
    boost::mt19937 gen_;
#endif
    REAL_DIST<double> dis_;
    NORM_DIST<double> insertion_length_gaussian_;
    INT_DIST<int> base_;
    INT_DIST<am_qual_t> quality_less_than_10_;
    INT_DIST<int> pos_on_read_1_;
    INT_DIST<int> pos_on_read_2_;
    INT_DIST<int> pos_on_read_1_not_head_and_tail_;
    INT_DIST<int> pos_on_read_2_not_head_and_tail_;
#elif defined(USE_ONEMKL_RANDOM)
    VSLStreamStatePtr stream_;
    double pe_frag_dist_mean_;
    double pe_frag_dist_std_dev_;
#endif
    am_read_len_t read_len_1_;
    am_read_len_t read_len_2_;
    am_read_len_t read_len_max_ { 0 };
};
} // namespace labw::art_modern
