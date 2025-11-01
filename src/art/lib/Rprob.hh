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

#include "art_modern_config.h"

#include "libam_support/utils/class_macros_utils.hh"

#if defined(USE_STL_RANDOM)
#include <random>
#elif defined(USE_BOOST_RANDOM)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#elif defined(USE_ONEMKL_RANDOM)
#include <mkl.h> // NOLINT
#elif defined(USE_PCG_RANDOM)
#include "libam_support/ds/pcg_32_c.hh"

#include <random>
#else
#error "Define USE_STL_RANDOM, USE_BOOST_RANDOM, USE_ONEMKL_RANDOM, USE_PCG_RANDOM for random generators!"
#endif

#include <cstddef>
#include <cstdint>
#include <vector>

namespace labw::art_modern {

class Rprob {
public:
    Rprob(double pe_frag_dist_mean, double pe_frag_dist_std_dev, int read_length);
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
     * Populate tmp_probs_ with read_len random values using r_prob().
     */
    void r_probs();
    /**
     * Generate an insertion length based on Gaussian distribution.
     * @return
     */
    int insertion_length();
    /**
     * Generate one of A, C, G, T.
     */
    char rand_base();
    /**
     * Refill tmp_qual_dists_.
     */
    void rand_quality_dist();
    int rand_quality_less_than_10();
    int rand_pos_on_read();
    int rand_pos_on_read_not_head_and_tail();
    int randint(int min, int max);
    std::vector<double> tmp_probs_;
    std::vector<int> tmp_qual_dists_;

private:
    static std::uint64_t seed();
    void public_init_();

#if defined(USE_ONEMKL_RANDOM)
    std::vector<int> cached_rand_pos_on_read_;
    std::size_t cached_rand_pos_on_read_index_ = 0;

    std::vector<int> cached_rand_pos_on_read_not_head_and_tail_;
    std::size_t cached_rand_pos_on_read_not_head_and_tail_index_ = 0;

    std::vector<double> cached_insertion_lengths_;
    std::size_t cached_insertion_lengths_index_ = 0;

    std::vector<int> cached_rand_base_indices_;
    std::size_t cached_rand_base_indices_index_ = 0;

    std::vector<int> cached_rand_quality_less_than_10_;
    std::size_t cached_rand_quality_less_than_10_index_ = 0;

    constexpr static std::size_t CACHE_SIZE_ = 4096;
#endif

#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM)
#if defined(USE_STL_RANDOM)
    std::mt19937 gen_;
#else
    pcg32_c gen_;
#endif
    std::uniform_real_distribution<double> dis_;
    std::normal_distribution<double> insertion_length_gaussian_;
    std::uniform_int_distribution<int> base_;
    std::uniform_int_distribution<int> strand_;
    std::uniform_int_distribution<int> quality_less_than_10_;
    std::uniform_int_distribution<int> quality_;
    std::uniform_int_distribution<int> pos_on_read_;
    std::uniform_int_distribution<int> pos_on_read_not_head_and_tail_;
#elif defined(USE_BOOST_RANDOM)
    boost::mt19937 gen_;
    boost::random::uniform_01<double> dis_;
    boost::random::normal_distribution<double> insertion_length_gaussian_;
    boost::random::uniform_int_distribution<int> base_;
    boost::random::uniform_int_distribution<int> strand_;
    boost::random::uniform_int_distribution<int> quality_less_than_10_;
    boost::random::uniform_int_distribution<int> quality_;
    boost::random::uniform_int_distribution<int> pos_on_read_;
    boost::random::uniform_int_distribution<int> pos_on_read_not_head_and_tail_;
#elif defined(USE_ONEMKL_RANDOM)
    VSLStreamStatePtr stream_;
    double pe_frag_dist_mean_;
    double pe_frag_dist_std_dev_;
#endif
    int read_length_;
};
} // namespace labw::art_modern
