#include "random_generator.hh"
#include "ArtConstants.hh"

namespace labw {
namespace art_modern {
#if defined(USE_STD_RANDOM)
    Rprob::Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev)
        : gen_(rd_())
        , dis_(0.0f, 1.0f)
        , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
        , base_(0, 3)
        , strand_(0, 1)
        , quality_less_than_10_(0, 10)
        , quality_(0, MAX_DIST_NUMBER)
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
}
#elif defined(USE_BOOST_RANDOM)
    Rprob::Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev)
        : gen_(rd_())
        , dis_(0.0f, 1.0f)
        , insertion_length_gaussian_(pe_frag_dist_mean, pe_frag_dist_std_dev)
        , base_(0, 3)
        , strand_(0, 1)
        , quality_less_than_10_(0, 10)
        , quality_(0, MAX_DIST_NUMBER)
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
}
#elif defined(USE_ONEMKL_RANDOM)
// Rprob::Rprob(float pe_frag_dist_mean, float pe_frag_dist_std_dev)
//     : rd_(), queue_(), gen_(queue_), dis_(0.0f, 1.0f)
//     , insertion_length_gaussian_(
//         pe_frag_dist_mean, pe_frag_dist_std_dev)
//     , base_(0, 3)
//     , strand_(0, 1)
//     , quality_less_than_10_(0, 10)
//     , quality_(0, MAX_DIST_NUMBER)
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
// }
#endif
}