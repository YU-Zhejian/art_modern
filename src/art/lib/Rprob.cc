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

#include "art/lib/Rprob.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/arithmetic_utils.hh"
#include "libam_support/utils/rand_utils.hh"

#if defined(USE_STL_RNDOM)
#include <random>
#endif

#if defined(USE_BOOST_RANDOM)
#include <boost/random/uniform_int_distribution.hpp>
#endif

#if defined(USE_ONEMKL_RANDOM)
#include <mkl.h>
#endif

#include <algorithm> // NOLINT for std::generate_n
#include <cstddef>
#include <vector>

namespace labw::art_modern {

void Rprob::public_init_()
{
    read_len_max_ = am_max(read_len_1_, read_len_2_);
    tmp_probs_.resize(read_len_max_);
}

#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM) || defined(USE_BOOST_RANDOM)
Rprob::~Rprob() = default;
Rprob::Rprob(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const am_read_len_t read_len_1,
    am_read_len_t read_len_2)
    : gen_(rand_seed())
    , dis_(0.0, 1.0)
    , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
    , base_(0, 3)
    , quality_less_than_10_(1, 10)
    , pos_on_read_1_(0, read_len_1 - 1)
    , pos_on_read_2_(0, read_len_2 - 1)
    , pos_on_read_1_not_head_and_tail_(1, read_len_1 - 2)
    , pos_on_read_2_not_head_and_tail_(1, read_len_2 - 2)
    , read_len_1_(read_len_1)
    , read_len_2_(read_len_2)
{
    public_init_();
}
#elif defined(USE_ONEMKL_RANDOM)
Rprob::~Rprob() { vslDeleteStream(&stream_); }

Rprob::Rprob(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const am_read_len_t read_len_1,
    const am_read_len_t read_len_2)
    : stream_() // Initialized in the function body
    , pe_frag_dist_mean_(pe_frag_dist_mean)
    , pe_frag_dist_std_dev_(pe_frag_dist_std_dev)
    , read_len_1_(read_len_1)
    , read_len_2_(read_len_2)
{
    vslNewStream(&stream_, VSL_BRNG_SFMT19937, rand_seed());
    public_init_();
}
#endif

void Rprob::r_probs(const std::size_t n)
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    std::generate_n(tmp_probs_.begin(), n, [this]() { return dis_(gen_); });
#elif defined(USE_ONEMKL_RANDOM)
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, static_cast<MKL_INT>(n), tmp_probs_.data(), 0.0, 1.0);
#endif
}

int Rprob::insertion_length()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return static_cast<int>(insertion_length_gaussian_(gen_));
#elif defined(USE_ONEMKL_RANDOM)
    if (cached_insertion_lengths_index_ == 0) {
        vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream_, CACHE_SIZE_, cached_insertion_lengths_.data(),
            pe_frag_dist_mean_, pe_frag_dist_std_dev_);
        cached_insertion_lengths_index_ = CACHE_SIZE_;
    }
    return static_cast<int>(cached_insertion_lengths_[--cached_insertion_lengths_index_]);
#endif
}

char Rprob::rand_base()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return ART_ACGT[base_(gen_)];
#elif defined(USE_ONEMKL_RANDOM)
    if (cached_rand_base_indices_index_ == 0) {
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, CACHE_SIZE_, cached_rand_base_indices_.data(), 0, 4);
        cached_rand_base_indices_index_ = CACHE_SIZE_;
    }
    return ART_ACGT[cached_rand_base_indices_[--cached_rand_base_indices_index_]];
#endif
}

am_qual_t Rprob::rand_quality_less_than_10()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return quality_less_than_10_(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    if (cached_rand_quality_less_than_10_index_ == 0) {
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, CACHE_SIZE_,
            reinterpret_cast<int*>(cached_rand_quality_less_than_10_.data()), 1, 10);
        cached_rand_quality_less_than_10_index_ = CACHE_SIZE_;
    }
    return cached_rand_quality_less_than_10_[--cached_rand_quality_less_than_10_index_];
#endif
}

int Rprob::rand_pos_on_read(const bool is_read1)
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return (is_read1 ? pos_on_read_1_ : pos_on_read_2_)(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    auto cached_rand_pos_on_read_index_
        = is_read1 ? cached_rand_pos_on_read_1_index_ : cached_rand_pos_on_read_2_index_;
    auto& cached_rand_pos_on_read_ = is_read1 ? cached_rand_pos_on_read_1_ : cached_rand_pos_on_read_2_;
    const auto read_length_ = is_read1 ? read_len_1_ : read_len_2_;
    if (cached_rand_pos_on_read_index_ == 0) {
        viRngUniform(
            VSL_RNG_METHOD_UNIFORM_STD, stream_, CACHE_SIZE_, cached_rand_pos_on_read_.data(), 0, read_length_);
        cached_rand_pos_on_read_index_ = CACHE_SIZE_;
    }
    return cached_rand_pos_on_read_[--cached_rand_pos_on_read_index_];
#endif
}
int Rprob::rand_pos_on_read_not_head_and_tail(const bool is_read1)
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return (is_read1 ? pos_on_read_1_not_head_and_tail_ : pos_on_read_2_not_head_and_tail_)(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    auto cached_rand_pos_on_read_not_head_and_tail_index_ = is_read1
        ? cached_rand_pos_on_read_1_not_head_and_tail_index_
        : cached_rand_pos_on_read_2_not_head_and_tail_index_;
    auto& cached_rand_pos_on_read_not_head_and_tail_
        = is_read1 ? cached_rand_pos_on_read_1_not_head_and_tail_ : cached_rand_pos_on_read_2_not_head_and_tail_;
    const auto read_length_ = is_read1 ? read_len_1_ : read_len_2_;
    if (cached_rand_pos_on_read_not_head_and_tail_index_ == 0) {
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, CACHE_SIZE_,
            cached_rand_pos_on_read_not_head_and_tail_.data(), 1, read_length_ - 1);
        cached_rand_pos_on_read_not_head_and_tail_index_ = CACHE_SIZE_;
    }
    return cached_rand_pos_on_read_not_head_and_tail_[--cached_rand_pos_on_read_not_head_and_tail_index_];
#endif
}
int Rprob::randint(const int min, const int max)
{
#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM)
    return std::uniform_int_distribution<int>(min, max - 1)(gen_);
#elif defined(USE_BOOST_RANDOM)
    return boost::random::uniform_int_distribution<int>(min, max - 1)(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    int result = 0;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, min, max);
    return result;
#endif
}

} // namespace labw::art_modern
