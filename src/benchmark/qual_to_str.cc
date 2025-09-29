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
 **/

#include "libam_support/Dtypes.hh"
#include "libam_support/utils/seq_utils.hh"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

using namespace labw::art_modern; // NOLINT

namespace {
void bench(const int run_times, const int rlen)
{
    std::cout << "run_times: " << run_times << " rlen: " << rlen << std::endl;
    std::vector<am_qual_t> q;
    q.reserve(rlen);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution dis(0, 40);
    for (int i = 0; i < rlen; i++) {
        q.emplace_back(dis(gen));
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_avx2(q.data(), rlen); // NOLINT
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "AVX2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_sse2(q.data(), rlen); // NOLINT
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "SSE2: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_mmx(q.data(), rlen); // NOLINT
    }

    end = std::chrono::high_resolution_clock::now();
    std::cout << "MMX: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str(q.data(), rlen); // NOLINT
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Scala (for loop): " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < run_times; i++) {
        volatile auto s = qual_to_str_foreach(q.data(), rlen); // NOLINT
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Scala (std::for_each): " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << std::endl;

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
