#pragma once

#include "art_modern_config.h"

#include "libam/utils/class_macros_utils.hh"

#if defined(USE_STL_RANDOM)
#include <random>
#elif defined(USE_BOOST_RANDOM)
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#elif defined(USE_ONEMKL_RANDOM)
#include <mkl.h> // NOLINT
#elif defined(USE_GSL_RANDOM)
#include <gsl/gsl_rng.h>
#elif defined(USE_PCG_RANDOM)
#include <pcg_random.hpp>

#include <random>
#else
#error                                                                                                                 \
    "Define USE_STL_RANDOM, USE_BOOST_RANDOM, USE_ONEMKL_RANDOM, USE_PCG_RANDOM or USE_GSL_RANDOM for random generators!"
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
     * Generate uniform [0, 1)
     */
    double r_prob();
    /**
     * Populate tmp_probs_ with n random values using r_prob().
     */
    void r_probs(std::size_t n);
    /**
     * Populate tmp_probs_ with read_len random values using r_prob().
     */
    void r_probs();
    int insertion_length();
    /**
     * Generate one of A, C, G, T.
     */
    char rand_base();
    void rand_quality_dist();
    int rand_quality_less_than_10();
    int rand_pos_on_read();
    int rand_pos_on_read_not_head_and_tail();
    int randint(int min, int max);
    std::vector<double> tmp_probs_;
    std::vector<int> tmp_qual_dists_;

private:
    static uint64_t seed();
    void public_init_();
#if defined(USE_STL_RANDOM) || defined(USE_PCG_RANDOM)
#if defined(USE_STL_RANDOM)
    std::mt19937 gen_;
#else
    pcg32_fast gen_;
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
#elif defined(USE_GSL_RANDOM)
    gsl_rng* r;
    double pe_frag_dist_mean_;
    double pe_frag_dist_std_dev_;
#endif
    int read_length_;
};
} // namespace labw::art_modern
