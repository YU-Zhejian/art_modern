#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include <mkl.h>

#include <gsl/gsl_rng.h>

#include <boost/random.hpp>

#include <absl/random/random.h>
#include <absl/random/uniform_int_distribution.h>
#include <absl/random/uniform_real_distribution.h>

#include "rprobs.hh"

class SlimRprobs {
public:
    virtual std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) = 0;
    virtual std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) = 0;
    virtual ~SlimRprobs() = default;
    virtual std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) = 0;
};

/**
 * Very, very, very slow.
 */
class [[maybe_unused]] SlimRprobsTrng : public SlimRprobs {
public:
    SlimRprobsTrng() = default;
    ~SlimRprobsTrng() override = default;
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
    std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gen_(); });
        return tmp_qual_dists_;
    }

private:
    std::random_device gen_;
};

class SlimRprobsStdRand : public SlimRprobs {
public:
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

    std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gen_(); });
        return tmp_qual_dists_;
    }

private:
    std::mt19937 gen_;
};

class SlimRprobsMKL : public SlimRprobs {
public:
    SlimRprobsMKL() { vslNewStream(&vsl_stream_, VSL_BRNG_MT19937, seed()); }
    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, vsl_stream_, N_BASES, tmp_qual_dists_.data(), a, b);
        return tmp_qual_dists_;
    }
    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, vsl_stream_, N_BASES, tmp_qual_dists_.data(), a, b);
        return tmp_qual_dists_;
    }

    std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) override
    {
        viRngUniformBits(VSL_RNG_METHOD_UNIFORM_STD, vsl_stream_, N_BASES, tmp_qual_dists_.data());
        return tmp_qual_dists_;
    }

    ~SlimRprobsMKL() override { vslDeleteStream(&vsl_stream_); }

private:
    VSLStreamStatePtr vsl_stream_ = nullptr;
};

class SlimRprobsBoost : public SlimRprobs {
public:
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

    std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gen_(); });
        return tmp_qual_dists_;
    }

private:
    boost::mt19937 gen_;
};

class SlimRprobsAbsl : public SlimRprobs {
public:
    SlimRprobsAbsl() = default;
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

    std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gen_(); });
        return tmp_qual_dists_;
    }

private:
    absl::BitGen gen_;
};

class SlimRprobsGsl : public SlimRprobs {
public:
    SlimRprobsGsl()
        : r(gsl_rng_alloc(gsl_rng_mt19937))
    {
        gsl_rng_set(r, seed());
    }

    std::vector<double> gen_doubles(std::vector<double>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gsl_rng_uniform(r); });
        return tmp_qual_dists_;
    }
    std::vector<int> gen_ints(std::vector<int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gsl_rng_uniform_int(r, b); });
        return tmp_qual_dists_;
    }

    std::vector<unsigned int> gen_bits(std::vector<unsigned int>& tmp_qual_dists_) override
    {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [this]() { return gsl_rng_get(r); });
        return tmp_qual_dists_;
    }
    ~SlimRprobsGsl() override { gsl_rng_free(r); }

private:
    gsl_rng* r;
};

namespace {
void bench(std::unique_ptr<SlimRprobs> rprobs, const std::string& name)
{
    std::vector<double> tmp_qual_dists_double_(N_BASES);
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    rprobs->gen_doubles(tmp_qual_dists_double_);
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << name << ": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;

    std::vector<int> tmp_qual_dists_int_(N_BASES);
    start = std::chrono::high_resolution_clock::now();
    rprobs->gen_ints(tmp_qual_dists_int_);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Int: " << name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << std::endl;

    std::vector<unsigned int> tmp_qual_dists_bits_(N_BASES);
    start = std::chrono::high_resolution_clock::now();
    rprobs->gen_bits(tmp_qual_dists_bits_);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Bits: " << name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << std::endl;
}

} // namespace

int main()
{
    bench(std::make_unique<SlimRprobsStdRand>(), "std::random");
    bench(std::make_unique<SlimRprobsMKL>(), "MKL");
    bench(std::make_unique<SlimRprobsBoost>(), "Boost");
    bench(std::make_unique<SlimRprobsGsl>(), "GSL");
    bench(std::make_unique<SlimRprobsAbsl>(), "Absl");
    // bench(std::make_unique<SlimRprobsTrng>(), "Trng");
}
