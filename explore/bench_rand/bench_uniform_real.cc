#include "benchmark_utils.hh"
#include "rprobs.hh"

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <absl/random/uniform_real_distribution.h>

#include <pcg_random.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <random>
#include <vector>

/**
 * ON ARM MACHINE:
 * Double: boost::random::uniform_01:                gmean:      7,315; mean/sd:       7,317/182us
 * Float: boost::random::uniform_01:                 gmean:      9,516; mean/sd:       9,518/193us
 * Double: boost::random::uniform_real_distribution: gmean:     21,279; mean/sd:      21,281/272us
 * Float: boost::random::uniform_real_distribution:  gmean:     29,092; mean/sd:      29,094/402us
 * Double: absl::uniform_real_distribution:          gmean:     13,394; mean/sd:      13,394/143us
 * Float: absl::uniform_real_distribution:           gmean:     15,073; mean/sd:      15,075/229us
 * STL's too slow to finish
 *
 * On x86_64 MACHINE:
 *  Double:                  boost::random::uniform_01: gmean:      1,664; mean/sd:       1,672/175 us
 *  _Float:                  boost::random::uniform_01: gmean:      1,768; mean/sd:       1,773/138 us
 *  Double:   boost::random::uniform_real_distribution: gmean:      1,630; mean/sd:        1,632/79 us
 *  _Float:   boost::random::uniform_real_distribution: gmean:      1,685; mean/sd:        1,686/44 us
 *  Double:            absl::uniform_real_distribution: gmean:      2,517; mean/sd:        2,517/34 us
 *  _Float:            absl::uniform_real_distribution: gmean:      2,614; mean/sd:        2,615/69 us
 *  Double:             std::uniform_real_distribution: gmean:      4,435; mean/sd:        4,435/14 us
 *  _Float:             std::uniform_real_distribution: gmean:      2,772; mean/sd:        2,772/11 us
 */
struct DistInfo {
    const char* name;
    std::function<double(pcg32_fast&)> gen_double;
    std::function<float(pcg32_fast&)> gen_float;
};

int main()
{
    pcg32_fast engine { };

    std::vector<double> tmp_qual_dists_ { };
    tmp_qual_dists_.reserve(N_BASES);

    std::vector<std::size_t> times { };

    std::vector<DistInfo> const dists
        = { { "boost::random::uniform_01", [](pcg32_fast& e) { return boost::random::uniform_01<double> { }(e); },
                [](pcg32_fast& e) { return boost::random::uniform_01<float> { }(e); } },
              { "boost::random::uniform_real_distribution",
                  [](pcg32_fast& e) { return boost::random::uniform_real_distribution<double> { 0.0, 1.0 }(e); },
                  [](pcg32_fast& e) { return boost::random::uniform_real_distribution<float> { 0.0, 1.0 }(e); } },
              { "absl::uniform_real_distribution",
                  [](pcg32_fast& e) { return absl::uniform_real_distribution<double> { 0.0, 1.0 }(e); },
                  [](pcg32_fast& e) { return absl::uniform_real_distribution<float> { 0.0, 1.0 }(e); } },
              { "std::uniform_real_distribution",
                  [](pcg32_fast& e) { return std::uniform_real_distribution<double> { 0.0, 1.0 }(e); },
                  [](pcg32_fast& e) { return std::uniform_real_distribution<float> { 0.0, 1.0 }(e); } } };

    for (const auto& dist : dists) {
        for (auto type : { std::make_pair("Double", true), std::make_pair("_Float", false) }) {
            times.clear();
            for (std::size_t i = 0; i < N_REPLICA; ++i) {
                auto start = std::chrono::high_resolution_clock::now();
                if (type.second) {
                    for (std::size_t j = 0; j < N_TIMES; ++j) {
                        std::generate_n(
                            tmp_qual_dists_.begin(), N_BASES, [&engine, &dist]() { return dist.gen_double(engine); });
                    }
                } else {
                    for (std::size_t j = 0; j < N_TIMES; ++j) {
                        std::generate_n(
                            tmp_qual_dists_.begin(), N_BASES, [&engine, &dist]() { return dist.gen_float(engine); });
                    }
                }
                auto end = std::chrono::high_resolution_clock::now();
                times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
            }
            std::cout << std::setw(10) << type.first << ": " << std::setw(42) << dist.name << ": " << describe(times) << " us"
                      << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
