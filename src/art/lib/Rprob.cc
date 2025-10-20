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

#include "art/lib/ArtConstants.hh"

#include "libam_support/Constants.hh"
#include "libam_support/utils/rand_utils.hh"

#if defined(USE_GSL_RANDOM)
#include <gsl/gsl_randist.h>
#endif

#if defined(USE_STL_RNDOM)
#include <random>
#endif

#if defined(USE_BOOST_RANDOM)
#include <boost/random/uniform_int_distribution.hpp>
#endif

#if defined(USE_ONEMKL_RANDOM)
#include <mkl.h>
#endif

#include <algorithm> // NOLINT
#include <cstddef>
#include <cstdint>
#include <vector>

namespace labw::art_modern {

std::uint64_t Rprob::seed()
{
    return rand_seed();
}

void Rprob::public_init_()
{
    tmp_qual_dists_.resize(read_length_);
    tmp_probs_.resize(read_length_);
}
void Rprob::r_probs() { r_probs(read_length_); }

#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM)
Rprob::~Rprob() = default;
Rprob::Rprob(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const int read_length)
    : gen_(seed())
    , dis_(0.0, 1.0)
    , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
    , base_(0, 3)
    , strand_(0, 1)
    , quality_less_than_10_(1, 10)
    , quality_(1, MAX_DIST_NUMBER)
    , pos_on_read_(0, read_length - 1)
    , pos_on_read_not_head_and_tail_(1, read_length - 2)
    , read_length_(read_length)
{
    public_init_();
}
#elif defined(USE_BOOST_RANDOM)
Rprob::~Rprob() = default;
Rprob::Rprob(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const int read_length)
    : gen_(seed())
    , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
    , base_(0, 3)
    , strand_(0, 1)
    , quality_less_than_10_(1, 10)
    , quality_(1, MAX_DIST_NUMBER)
    , pos_on_read_(0, read_length - 1)
    , pos_on_read_not_head_and_tail_(1, read_length - 2)
    , read_length_(read_length)
{
    public_init_();
}
#elif defined(USE_ONEMKL_RANDOM)
Rprob::~Rprob() = default;
Rprob::Rprob(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const int read_length)
    : stream_() // Initialized in the function body
    , pe_frag_dist_mean_(pe_frag_dist_mean)
    , pe_frag_dist_std_dev_(pe_frag_dist_std_dev)
    , read_length_(read_length)
{
    vslNewStream(&stream_, VSL_BRNG_MT19937, seed());
    public_init_();
}
#elif defined(USE_GSL_RANDOM)
Rprob::Rprob(double pe_frag_dist_mean, double pe_frag_dist_std_dev, int read_length)
    : pe_frag_dist_mean_(pe_frag_dist_mean)
    , pe_frag_dist_std_dev_(pe_frag_dist_std_dev)
    , read_length_(read_length)
{
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed());
    public_init_();
}
Rprob::~Rprob() { gsl_rng_free(r); }
#endif

double Rprob::r_prob()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return dis_(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    double result = 0.0;
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 0.0, 1.0);
    return result;
#elif defined(USE_GSL_RANDOM)
    return gsl_rng_uniform(r);
#endif
}

void Rprob::r_probs(const std::size_t n)
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_GSL_RANDOM) || defined(USE_PCG_RANDOM)
    std::generate_n(tmp_probs_.begin(), n, [this]() { return r_prob(); });
#elif defined(USE_ONEMKL_RANDOM)
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, static_cast<MKL_INT>(n), tmp_probs_.data(), 0.0, 1.0);
#endif
}

int Rprob::insertion_length()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return static_cast<int>(insertion_length_gaussian_(gen_));
#elif defined(USE_ONEMKL_RANDOM)
    double result = 0.0;
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream_, 1, &result, pe_frag_dist_mean_, pe_frag_dist_std_dev_);
    return static_cast<int>(result);
#elif defined(USE_GSL_RANDOM)
    return static_cast<int>(gsl_ran_gaussian(r, pe_frag_dist_std_dev_) + pe_frag_dist_mean_);
#endif
}

char Rprob::rand_base()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return ART_ACGT[base_(gen_)];
#elif defined(USE_ONEMKL_RANDOM)
    int index = 0;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &index, 0, 4);
    return ART_ACGT[index];
#elif defined(USE_GSL_RANDOM)
    return ART_ACGT[gsl_rng_uniform_int(r, 4)];
#endif
}

void Rprob::rand_quality_dist()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    std::generate_n(tmp_qual_dists_.begin(), read_length_, [this]() { return quality_(gen_); });
#elif defined(USE_ONEMKL_RANDOM)
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, read_length_, tmp_qual_dists_.data(), 1, MAX_DIST_NUMBER);
#elif defined(USE_GSL_RANDOM)
    std::generate_n(tmp_qual_dists_.begin(), read_length_,
        [this]() { return static_cast<int>(gsl_rng_uniform_int(r, MAX_DIST_NUMBER - 1) + 1); });
#endif
}

int Rprob::rand_quality_less_than_10()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return quality_less_than_10_(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    int result = 0;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 1, 10);
    return result;
#elif defined(USE_GSL_RANDOM)
    return static_cast<int>(gsl_rng_uniform_int(r, 9) + 1);
#endif
}

int Rprob::rand_pos_on_read()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return pos_on_read_(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    int result = 0;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 0, read_length_);
    return result;
#elif defined(USE_GSL_RANDOM)
    return static_cast<int>(gsl_rng_uniform_int(r, read_length_));
#endif
}
int Rprob::rand_pos_on_read_not_head_and_tail()
{
#if defined(USE_STL_RANDOM) || defined(USE_BOOST_RANDOM) || defined(USE_PCG_RANDOM)
    return pos_on_read_not_head_and_tail_(gen_);
#elif defined(USE_ONEMKL_RANDOM)
    int result = 0;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 1, read_length_ - 1);
    return result;
#elif defined(USE_GSL_RANDOM)
    return static_cast<int>(gsl_rng_uniform_int(r, read_length_ - 2) + 1);
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
#elif defined(USE_GSL_RANDOM)
    return static_cast<int>(gsl_rng_uniform_int(r, max - min) + min);
#endif
}

} // namespace labw::art_modern
