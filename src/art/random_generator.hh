#pragma once

#include "art_modern_config.h"

#include "libam/utils/class_macros_utils.hh"

#if defined(USE_STL_RANDOM)
#include <random>
#elif defined(USE_BOOST_RANDOM)
#include <boost/random.hpp>
#elif defined(USE_ONEMKL_RANDOM)
#include <mkl.h> // NOLINT
#elif defined(USE_GSL_RANDOM)
#include <gsl/gsl_rng.h>
#else
#error "Define USE_STL_RANDOM, USE_BOOST_RANDOM, USE_ONEMKL_RANDOM or USE_GSL_RANDOM for random generators!"
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

    double r_prob();
    void r_probs(std::size_t n);
    void r_probs();
    int insertion_length();
    char rand_base();
    void rand_quality();
    int rand_quality_less_than_10();
    int rand_pos_on_read();
    int rand_pos_on_read_not_head_and_tail();
    int randint(int min, int max);
    std::vector<double> tmp_probs_;
    std::vector<int> tmp_qual_dists_;

private:
    static uint64_t seed();
    void public_init_();
#if defined(USE_STL_RANDOM)
    std::mt19937 gen_;
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
    boost::uniform_real<double> dis_;
    boost::normal_distribution<double> insertion_length_gaussian_;
    boost::uniform_smallint<int> base_;
    boost::uniform_smallint<int> strand_;
    boost::uniform_smallint<int> quality_less_than_10_;
    boost::uniform_int<int> quality_;
    boost::uniform_int<int> pos_on_read_;
    boost::uniform_int<int> pos_on_read_not_head_and_tail_;
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
