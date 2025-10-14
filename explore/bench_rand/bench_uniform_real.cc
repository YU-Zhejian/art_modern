#include "rprobs.hh"
#include "benchmark_utils.hh"

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <absl/random/uniform_real_distribution.h>

#include <pcg_random.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
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
 */

int main()
{
    pcg32_fast engine { };

    std::vector<double> tmp_qual_dists_ { };
    tmp_qual_dists_.reserve(N_BASES);

    std::vector<std::size_t> times { };

    boost::random::uniform_01<double> boost_dist { };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &boost_dist]() { return boost_dist(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Double: " << "boost::random::uniform_01" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    boost::random::uniform_01<float> boost_dist_f { };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &boost_dist_f]() { return boost_dist_f(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Float: " << "boost::random::uniform_01" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    boost::random::uniform_real_distribution<double> boost_dist01 { 0.0, 1.0 };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &boost_dist01]() { return boost_dist01(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Double: " << "boost::random::uniform_real_distribution" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();
    
    boost::random::uniform_real_distribution<float> boost_dist01_f { 0.0, 1.0 };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &boost_dist01_f]() { return boost_dist01_f(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Float: " << "boost::random::uniform_real_distribution" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    absl::uniform_real_distribution<double> absl_dist { 0.0, 1.0 };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &absl_dist]() { return absl_dist(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Double: " << "absl::uniform_real_distribution" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    absl::uniform_real_distribution<float> absl_dist_f { 0.0, 1.0 };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &absl_dist_f]() { return absl_dist_f(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Float: " << "absl::uniform_real_distribution" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    std::uniform_real_distribution<double> std_dist { 0.0, 1.0 };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    for(std::size_t j = 0; j < N_TIMES; ++j){
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &std_dist]() { return std_dist(engine); });}}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Double: " << "std::uniform_real_distribution" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    std::uniform_real_distribution<float> std_dist_f { 0.0, 1.0 };
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
            for(std::size_t j = 0; j < N_TIMES; ++j){std::generate_n(tmp_qual_dists_.begin(), N_BASES, [&engine, &std_dist_f]() { return std_dist_f(engine); });}
    auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << "Float: " << "std::uniform_real_distribution" << ": "
              << describe(times) << "us"
              << std::endl;
    times.clear();

    return EXIT_SUCCESS;
}
