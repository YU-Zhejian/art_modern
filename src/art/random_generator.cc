#include "random_generator.hh"
#include "ArtConstants.hh"

#include <algorithm>
#include <chrono>
#include <thread>
#include <vector>

#if defined(USE_GSL_RANDOM)
#include <gsl/gsl_randist.h>
#endif

namespace labw::art_modern {

long Rprob::seed()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
               .count()
        * static_cast<long>(std::hash<std::thread::id>()(std::this_thread::get_id()));
}

#if defined(USE_STL_RANDOM)
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
}

double Rprob::r_prob() { return dis_(gen_); }

int Rprob::insertion_length() { return static_cast<int>(insertion_length_gaussian_(gen_)); }

char Rprob::rand_base() { return ART_ACGT[base_(gen_)]; }

void Rprob::rand_quality(std::vector<int>& qual_dist)
{
    std::generate_n(qual_dist.begin(), read_length_, [this]() { return quality_(gen_); });
}

void Rprob::r_probs(std::vector<double>& result)
{
    std::generate_n(result.begin(), result.size(), [this]() { return r_prob(); });
}

int Rprob::rand_quality_less_than_10() { return quality_less_than_10_(gen_); }

int Rprob::rand_pos_on_read() { return pos_on_read_(gen_); }
int Rprob::rand_pos_on_read_not_head_and_tail() { return pos_on_read_not_head_and_tail_(gen_); }
int Rprob::randint(int min, int max) { return std::uniform_int_distribution<int>(min, max - 1)(gen_); }

#elif defined(USE_BOOST_RANDOM)
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
}

double Rprob::r_prob() { return dis_(gen_); }

void Rprob::r_probs(std::vector<double>& result)
{
    std::generate_n(result.begin(), result.size(), [this]() { return r_prob(); });
}

int Rprob::insertion_length() { return static_cast<int>(insertion_length_gaussian_(gen_)); }

char Rprob::rand_base() { return ART_ACGT[base_(gen_)]; }

void Rprob::rand_quality(std::vector<int>& qual_dist)
{
    std::generate_n(qual_dist.begin(), read_length_, [this]() { return quality_(gen_); });
}

int Rprob::rand_quality_less_than_10() { return quality_less_than_10_(gen_); }

int Rprob::rand_pos_on_read() { return pos_on_read_(gen_); }
int Rprob::rand_pos_on_read_not_head_and_tail() { return pos_on_read_not_head_and_tail_(gen_); }
int Rprob::randint(int min, int max) { return std::uniform_int_distribution<int>(min, max - 1)(gen_); }

#elif defined(USE_ONEMKL_RANDOM)
Rprob::~Rprob() = default;
Rprob::Rprob(const double pe_frag_dist_mean, const double pe_frag_dist_std_dev, const int read_length)
    : stream_() // Initialized in the function body
    , pe_frag_dist_mean_(pe_frag_dist_mean)
    , pe_frag_dist_std_dev_(pe_frag_dist_std_dev)
    , read_length_(read_length)
{
    vslNewStream(&stream_, VSL_BRNG_MT19937, seed());
}

double Rprob::r_prob()
{
    double result;
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 0.0, 1.0);
    return result;
}

void Rprob::r_probs(std::vector<double>& result)
{
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, result.size(), result.data(), 0.0, 1.0);
}

int Rprob::insertion_length()
{
    double result;
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream_, 1, &result, pe_frag_dist_mean_, pe_frag_dist_std_dev_);
    return static_cast<int>(result);
}

char Rprob::rand_base()
{
    unsigned int index;
    viRngUniformBits32(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &index);
    return "ACGT"[index & 0b11];
}

void Rprob::rand_quality(std::vector<int>& qual_dist)
{
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, read_length_, qual_dist.data(), 1, MAX_DIST_NUMBER);
}

int Rprob::rand_quality_less_than_10()
{
    int result;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 1, 10);
    return result;
}

int Rprob::rand_pos_on_read()
{
    int result;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 0, read_length_ - 1);
    return result;
}

int Rprob::rand_pos_on_read_not_head_and_tail()
{
    int result;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, 1, read_length_ - 2);
    return result;
}
int Rprob::randint(const int min, const int max)
{
    int result;
    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, 1, &result, min, max);
    return result;
}

#elif defined(USE_GSL_RANDOM)
Rprob::Rprob(double pe_frag_dist_mean, double pe_frag_dist_std_dev, int read_length)
    : pe_frag_dist_mean_(pe_frag_dist_mean)
    , pe_frag_dist_std_dev_(pe_frag_dist_std_dev)
    , read_length_(read_length)
{
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed());
}
Rprob::~Rprob() { gsl_rng_free(r); }

double Rprob::r_prob() { return gsl_rng_uniform(r); }

void Rprob::r_probs(std::vector<double>& result)
{
    std::generate_n(result.begin(), result.size(), [this]() { return r_prob(); });
}

int Rprob::insertion_length()
{
    return static_cast<int>(gsl_ran_gaussian(r, pe_frag_dist_std_dev_) + pe_frag_dist_mean_);
}

char Rprob::rand_base() { return ART_ACGT[gsl_rng_uniform_int(r, 4)]; }

void Rprob::rand_quality(std::vector<int>& qual_dist)
{
    std::generate_n(qual_dist.begin(), read_length_,
        [this]() { return static_cast<int>(gsl_rng_uniform_int(r, MAX_DIST_NUMBER) + 1); });
}

int Rprob::rand_quality_less_than_10() { return static_cast<int>(gsl_rng_uniform_int(r, 9) + 1); }

int Rprob::rand_pos_on_read() { return static_cast<int>(gsl_rng_uniform_int(r, read_length_)); }

int Rprob::rand_pos_on_read_not_head_and_tail()
{
    return static_cast<int>(gsl_rng_uniform_int(r, read_length_ - 2) + 1);
}

int Rprob::randint(int min, int max) { return static_cast<int>(gsl_rng_uniform_int(r, max - min) + min); }
#endif
}
