#include "random_generator.hh"
#include "ArtConstants.hh"

#if defined(USE_GSL_RANDOM)
#include <gsl/gsl_randist.h>
#endif

namespace labw {
namespace art_modern {
#if defined(USE_STL_RANDOM)
    Rprob::~Rprob() = default;
    Rprob::Rprob(const float pe_frag_dist_mean, const float pe_frag_dist_std_dev, const int read_length)
        : gen_()
        , dis_(0.0f, 1.0f)
        , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
        , base_(0, 3)
        , strand_(0, 1)
        , quality_less_than_10_(1, 10)
        , quality_(1, MAX_DIST_NUMBER)
        , pos_on_read_(0, read_length - 1)
        , pos_on_read_not_head_and_tail_(1, read_length - 2)
    {
    }

    double Rprob::r_prob() { return dis_(gen_); }

    int Rprob::insertion_length() { return static_cast<int>(insertion_length_gaussian_(gen_)); }

    char Rprob::rand_base()
    {
        switch (base_(gen_)) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        default:
            return 'T';
        }
    }

    int Rprob::rand_quality() { return quality_(gen_); }

    int Rprob::rand_quality_less_than_10() { return quality_less_than_10_(gen_); }

    int Rprob::rand_pos_on_read() { return pos_on_read_(gen_); }
    int Rprob::rand_pos_on_read_not_head_and_tail() { return pos_on_read_not_head_and_tail_(gen_); }

#elif defined(USE_BOOST_RANDOM)
    Rprob::~Rprob() = default;
    Rprob::Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev)
        : gen_()
        , dis_(0.0f, 1.0f)
        , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
        , base_(0, 3)
        , strand_(0, 1)
        , quality_less_than_10_(1, 10)
        , quality_(1, MAX_DIST_NUMBER)
    {
    }

    double Rprob::r_prob() { return dis_(gen_); }

    int Rprob::insertion_length() { return static_cast<int>(insertion_length_gaussian_(gen_)); }

    char Rprob::rand_base()
    {
        switch (base_(gen_)) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        default:
            return 'T';
        }
    }

    int Rprob::rand_quality() { return quality_(gen_); }

    int Rprob::rand_quality_less_than_10() { return quality_less_than_10_(gen_); }

#elif defined(USE_ONEMKL_RANDOM)
    Rprob::~Rprob() = default;
    // Rprob::Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev)
    //     : rd_(), queue_(), gen_(queue_), dis_(0.0f, 1.0f)
    //     , insertion_length_gaussian_(
    //         pe_frag_dist_mean, pe_frag_dist_std_dev)
    //     , base_(0, 3)
    //     , strand_(0, 1)
    //     , quality_less_than_10_(1, 10)
    //     , quality_(1, MAX_DIST_NUMBER)
    //{
    // }
    //
    // double Rprob::r_prob() { return dis_(gen_); }
    //
    // int Rprob::insertion_length() { return static_cast<int>(
    //     insertion_length_gaussian_(gen_));
    // }
    //
    // char Rprob::rand_base()
    //{
    //     int n = 1;
    //     sycl::usm_allocator<int, sycl::usm::alloc::shared> allocator(queue_);
    //     std::vector<int, decltype(allocator)> r(n, allocator);
    //     auto event = oneapi::mkl::rng::generate(base_, gen_, n, r.data());
    //     event.wait();
    //     switch (r[0]) {
    //         case 0:
    //             return 'A';
    //         case 1:
    //             return 'C';
    //         case 2:
    //             return 'G';
    //         default:
    //             return 'T';
    //     }
    // }
    //
    // int Rprob::rand_quality()
    //{
    //     return quality_(gen_);
    // }
    //
    // int Rprob::rand_quality_less_than_10()
    //{
    //     return quality_less_than_10_(gen_);
    // }

#elif defined(USE_GSL_RANDOM)

    Rprob::Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev, int read_length)
        : pe_frag_dist_mean_(pe_frag_dist_mean)
        , pe_frag_dist_std_dev_(pe_frag_dist_std_dev)
        , read_length_(read_length)
    {
        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    Rprob::~Rprob() { gsl_rng_free(r); }

    double Rprob::r_prob() { return gsl_rng_uniform(r); }

    int Rprob::insertion_length()
    {
        return static_cast<int>(gsl_ran_gaussian(r, pe_frag_dist_std_dev_) + pe_frag_dist_mean_);
    }

    char Rprob::rand_base()
    {
        switch (static_cast<int>(gsl_rng_uniform_int(r, 4))) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        default:
            return 'T';
        }
    }

    int Rprob::rand_quality() { return static_cast<int>(gsl_rng_uniform_int(r, MAX_DIST_NUMBER) + 1); }

    int Rprob::rand_quality_less_than_10() { return static_cast<int>(gsl_rng_uniform_int(r, 9) + 1); }

    int Rprob::rand_pos_on_read() { return gsl_rng_uniform_int(r, read_length_); }

    int Rprob::rand_pos_on_read_not_head_and_tail() { return gsl_rng_uniform_int(r, read_length_ - 2) + 1; }

#endif
}
}
