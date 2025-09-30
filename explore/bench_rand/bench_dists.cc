#include "class_utils.hh"
#include "gsl_rng_wrapper.hh"
#include "rprobs.hh"

#include <pcg_random.hpp>

#include <boost/random.hpp>

#include <absl/random/gaussian_distribution.h>
#include <absl/random/uniform_int_distribution.h>
#include <absl/random/uniform_real_distribution.h>

#include <gsl/gsl_randist.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <random>
#include <thread>
#include <utility>
#include <vector>

class SlimRprobs {
public:
    virtual std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) = 0;
    virtual std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) = 0;
    virtual std::vector<double> gen_normal(std::vector<double>& tmp_qual_dists_) = 0;
    SlimRprobs() = default;
    virtual ~SlimRprobs() = default;
    DELETE_COPY_MOVE(SlimRprobs);
};

class SlimRprobsStdRand : public SlimRprobs {
public:
    DELETE_COPY_MOVE(SlimRprobsStdRand);
    SlimRprobsStdRand()
        : gen_(seed())
    {
    }
    ~SlimRprobsStdRand() override = default;
    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        std::uniform_real_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        std::uniform_int_distribution<int> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<double> gen_normal(std::vector<double>& tmp_qual_dists_) override
    {
        std::normal_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

private:
    pcg32 gen_;
};

class SlimRprobsBoost : public SlimRprobs {
public:
    DELETE_COPY_MOVE(SlimRprobsBoost);
    SlimRprobsBoost()
        : gen_(seed())
    {
    }
    ~SlimRprobsBoost() override = default;
    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        boost::uniform_real<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        boost::uniform_int<int> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<double> gen_normal(std::vector<double>& tmp_qual_dists_) override
    {
        boost::normal_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

private:
    pcg32 gen_;
};

class SlimRprobsAbsl : public SlimRprobs {
public:
    DELETE_COPY_MOVE(SlimRprobsAbsl);
    SlimRprobsAbsl()
        : gen_(seed())
    {
    }
    ~SlimRprobsAbsl() override = default;
    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        absl::uniform_real_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        absl::uniform_int_distribution<int> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<double> gen_normal(std::vector<double>& tmp_qual_dists_) override
    {
        absl::gaussian_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

private:
    pcg32 gen_;
};

class SlimRprobsGslNative : public SlimRprobs {
public:
    DELETE_COPY_MOVE(SlimRprobsGslNative);
    SlimRprobsGslNative()
        : gen_(gsl_rng_alloc(gsl_rng_mt19937))
    {
        gsl_rng_set(gen_, seed());
    }
    ~SlimRprobsGslNative() override { gsl_rng_free(gen_); }
    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gsl_rng_uniform(gen_) * (b - a) + a; });
        return tmp_qual_dists_;
    }

    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gsl_rng_uniform_int(gen_, b + 1) + a; });
        return tmp_qual_dists_;
    }

    std::vector<double> gen_normal(std::vector<double>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gsl_ran_gaussian(gen_, b) + a; });
        return tmp_qual_dists_;
    }

private:
    gsl_rng* gen_;
};

class SlimRprobsGslBoost : public SlimRprobs {
public:
    DELETE_COPY_MOVE(SlimRprobsGslBoost);
    SlimRprobsGslBoost()
        : gen_(gsl_rng_mt19937, seed())
    {
    }
    ~SlimRprobsGslBoost() override = default;
    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        boost::random::uniform_real_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        boost::random::uniform_int_distribution<int> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

    std::vector<double> gen_normal(std::vector<double>& tmp_qual_dists_) override
    {
        boost::normal_distribution<double> dist(a, b);
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&dist, this]() { return dist(gen_); });
        return tmp_qual_dists_;
    }

private:
    GslRngWrapper gen_;
};

namespace {
void bench(std::unique_ptr<SlimRprobs> rprobs, const std::string& name)
{
    std::vector<double> tmp_qual_dists_double_(N_BASES);
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N_TIMES; i++) {
        rprobs->gen_doubles(tmp_qual_dists_double_);
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << name << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    std::vector<int> tmp_qual_dists_int_(N_BASES);
    start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N_TIMES; i++) {
        rprobs->gen_ints(tmp_qual_dists_int_);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Int: " << name << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    std::vector<double> tmp_qual_dists_normal_(N_BASES);
    start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N_TIMES; i++) {
        rprobs->gen_normal(tmp_qual_dists_normal_);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Normal: " << name << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;
}

} // namespace

int main()
{
    bench(std::make_unique<SlimRprobsStdRand>(), "std::random");
    bench(std::make_unique<SlimRprobsBoost>(), "Boost");
    bench(std::make_unique<SlimRprobsAbsl>(), "Absl");

    bench(std::make_unique<SlimRprobsGslNative>(), "GSL native");
    bench(std::make_unique<SlimRprobsGslBoost>(), "GSL Boost");
    return EXIT_SUCCESS;
}
