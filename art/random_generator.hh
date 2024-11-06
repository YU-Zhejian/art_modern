#pragma once
#include "art_modern_config.h"
#if defined(USE_STL_RANDOM)
#include <random>
#elif defined(USE_BOOST_RANDOM)
#include <boost/random.hpp>
#error "Under construction."
#elif defined(USE_ONEMKL_RANDOM)
#include <oneapi/mkl/rng.hpp>
#include <sycl/sycl.hpp>
#error "Under construction."
#elif defined(USE_GSL_RANDOM)
#include <gsl/gsl_rng.h>
#else
#error "Define USE_STL_RANDOM, USE_BOOST_RANDOM, USE_ONEMKL_RANDOM or USE_GSL_RANDOM for random generators!"
#endif

namespace labw {
namespace art_modern {

    class Rprob {
    public:
        Rprob(double pe_frag_dist_mean, double pe_frag_dist_std_dev, int read_length);
        double r_prob();
        int insertion_length();
        char rand_base();
        int rand_quality();
        int rand_quality_less_than_10();
        ~Rprob();
        int rand_pos_on_read();
        int rand_pos_on_read_not_head_and_tail();

    private:
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
        boost::uniform_int_distribution<int> pos_on_read_;
        boost::uniform_int_distribution<int> pos_on_read_not_head_and_tail_;
#elif defined(USE_ONEMKL_RANDOM)
#error "TODO"
//    sycl::device rd_;
//    sycl::queue queue_;
//    oneapi::mkl::rng::mt19937 gen_;
//    oneapi::mkl::rng::uniform<double, oneapi::mkl::rng::uniform_method::by_default> dis_;
//    oneapi::mkl::rng::gaussian<double, oneapi::mkl::rng::gaussian_method::by_default> insertion_length_gaussian_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> base_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> strand_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> quality_less_than_10_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> quality_;
#elif defined(USE_GSL_RANDOM)
        const gsl_rng_type* T;
        gsl_rng* r;
        double pe_frag_dist_mean_;
        double pe_frag_dist_std_dev_;
        int read_length_;
#endif
    };
}
}
