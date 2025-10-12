/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * TODO: Use Geometric mean to better represent the performance
 **/

#include "benchmark_utils.hh"

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/seq_utils.hh"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace labw::art_modern; // NOLINT

constexpr std::size_t N_REPLICA = 200UL;

namespace {

void bech_impl(std::string (*f)(const signed char*, std::size_t), const std::string& name, const int run_times,
    const int rlen, const std::vector<am_qual_t>& q)
{
    std::vector<std::size_t> times;
    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < run_times; i++) {
            volatile auto s = f(q.data(), rlen); // NOLINT
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << std::setw(25) << (name + ": ") << describe(times) << std::endl;
}

void bench(const int run_times, const int rlen)
{
    std::cout << ">run_times: " << run_times << " rlen: " << rlen << std::endl;
    std::vector<am_qual_t> q;
    q.reserve(rlen);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution dis(0, 40);
    for (int i = 0; i < rlen; i++) {
        q.emplace_back(dis(gen));
    }

    bech_impl(qual_to_str_for_loop, "Scala (for loop)", run_times, rlen, q);
    bech_impl(qual_to_str_foreach, "Scala (std::for_each)", run_times, rlen, q);
    bech_impl(qual_to_str_mmx, "MMX", run_times, rlen, q);
    bech_impl(qual_to_str_sse2, "SSE2", run_times, rlen, q);
    bech_impl(qual_to_str_avx2, "AVX2", run_times, rlen, q);

    if (qual_to_str_sse2(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        std::abort();
    }

    if (qual_to_str_mmx(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        std::abort();
    }

    if (qual_to_str_for_loop(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        std::abort();
    }

    if (qual_to_str_avx2(q.data(), rlen) != qual_to_str_foreach(q.data(), rlen)) {
        std::abort();
    }
}

} // namespace
int main()
{
    bench(1000000, 36);
    bench(1000000, 100);
    bench(1000000, 150);
    bench(1000000, 200);
    bench(1000000, 300);
    bench(100000, 1000);
    bench(10000, 10000);
    return 0;
}
