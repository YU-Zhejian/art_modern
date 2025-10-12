#include "rprobs.hh"

#include <boost/random/uniform_real_distribution.hpp>

#include <absl/random/uniform_real_distribution.h>

#include <pcg_random.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

namespace {
// FIXME: Who gives me this code?
template <typename Engine> double pcg32_generate_uniform_double(Engine& engine, const double min, const double max)
{
    return ldexpl((static_cast<uint64_t>(engine()) << 32) | engine(), -64) * (max - min) + min;
}

template <typename Engine> float pcg32_generate_uniform_float(Engine& engine, const float min, const float max)
{
    return ldexpf(engine(), -32) * (max - min) + min;
}

} // namespace
int main()
{
    pcg32_fast engine {};

    std::vector<double> tmp_qual_dists_ {};
    tmp_qual_dists_.reserve(N_BASES);
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;

    start = std::chrono::high_resolution_clock::now();
    boost::random::uniform_real_distribution<double> boost_dist { 0.0, 1.0 };
    for (int i = 0; i < N_TIMES; ++i) {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &boost_dist]() { return boost_dist(engine); });
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << "boost::random::uniform_real_distribution" << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    start = std::chrono::high_resolution_clock::now();
    absl::uniform_real_distribution<double> absl_dist { 0.0, 1.0 };
    for (int i = 0; i < N_TIMES; ++i) {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &absl_dist]() { return absl_dist(engine); });
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << "absl::uniform_real_distribution" << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::uniform_real_distribution<double> std_dist { 0.0, 1.0 };
    for (int i = 0; i < N_TIMES; ++i) {
        std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &std_dist]() { return std_dist(engine); });
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << "std::uniform_real_distribution" << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_TIMES; ++i) {
        std::generate_n(
            tmp_qual_dists_.begin(), N_BASES, [&engine]() { return pcg32_generate_uniform_double(engine, 0.0, 1.0); });
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << "pcg32_generate_uniform_double" << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_TIMES; ++i) {
        std::generate_n(
            tmp_qual_dists_.begin(), N_BASES, [&engine]() { return pcg32_generate_uniform_float(engine, 0.0, 1.0); });
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Double: " << "pcg32_generate_uniform_float" << ": "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) << "ns"
              << std::endl;

    return EXIT_SUCCESS;
}
