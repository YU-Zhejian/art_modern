#pragma once

#if defined(USE_STD_RANDOM)
#include <random>
#elif defined(USE_BOOST_RANDOM)
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#elif defined(USE_ONEMKL_RANDOM)
#include <oneapi/mkl/rng.hpp>
#include <sycl/sycl.hpp>
#elif defined(USE_GSL_RANDOM)
#include <gsl/gsl_rng.h>
#else
#error "Define USE_STD_RANDOM, USE_BOOST_RANDOM for random generators!"
#endif

namespace labw {
namespace art_modern {

    class Rprob {
    public:
        explicit Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev);
        double r_prob();
        int insertion_length();
        char rand_base();
        int rand_quality();
        int rand_quality_less_than_10();
        ~Rprob();

    private:
#if defined(USE_STD_RANDOM)
        std::mt19937 gen_;
        std::uniform_real_distribution<float> dis_;
        std::normal_distribution<float> insertion_length_gaussian_;
        std::uniform_int_distribution<int> base_;
        std::uniform_int_distribution<int> strand_;
        std::uniform_int_distribution<int> quality_less_than_10_;
        std::uniform_int_distribution<int> quality_;
#elif defined(USE_BOOST_RANDOM)
        boost::mt19937 gen_;
        boost::uniform_real<float> dis_;
        boost::normal_distribution<float> insertion_length_gaussian_;
        boost::uniform_smallint<int> base_;
        boost::uniform_smallint<int> strand_;
        boost::uniform_smallint<int> quality_less_than_10_;
        boost::uniform_int<int> quality_;
#elif defined(USE_ONEMKL_RANDOM)
#error "TODO"
//    sycl::device rd_;
//    sycl::queue queue_;
//    oneapi::mkl::rng::mt19937 gen_;
//    oneapi::mkl::rng::uniform<float, oneapi::mkl::rng::uniform_method::by_default> dis_;
//    oneapi::mkl::rng::gaussian<float, oneapi::mkl::rng::gaussian_method::by_default> insertion_length_gaussian_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> base_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> strand_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> quality_less_than_10_;
//    oneapi::mkl::rng::uniform<int, oneapi::mkl::rng::uniform_method::by_default> quality_;
#elif defined(USE_GSL_RANDOM)
        const gsl_rng_type* T;
        gsl_rng* r;
        float pe_frag_dist_mean_;
        float pe_frag_dist_std_dev_;
#endif
    };
}
}