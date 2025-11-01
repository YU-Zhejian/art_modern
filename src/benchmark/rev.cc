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

#include "benchmark_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/seq/qreverse.h"
#include "libam_support/utils/seq_utils.hh"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace labw::art_modern; // NOLINT

constexpr std::size_t N_REPLICA = 2000UL;
constexpr std::size_t BYTES_TO_PROCESS = 20ULL * M_SIZE;

namespace {

template <typename T>
void bech_impl(void (*rev_impl)(void*, std::size_t), const std::string& name, const std::size_t run_times,
    const int rlen, const std::vector<std::vector<T>>& qs)
{
    std::vector<std::vector<T>> qs_copy_impl;
    std::vector<std::vector<T>> qs_copy_correct;
    for (const auto& q : qs) {
        qs_copy_impl.emplace_back(q);
        qs_copy_correct.emplace_back(q);
    }
    std::vector<std::size_t> times_impl;
    std::vector<std::size_t> times_correct;

    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < run_times; i++) {
            rev_impl(static_cast<void*>(qs_copy_impl[j].data()), rlen);
        }
        auto end = std::chrono::high_resolution_clock::now();
        times_impl.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());

        start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < run_times; i++) {
            reverse(qs_copy_correct[j].data(), rlen);
        }
        end = std::chrono::high_resolution_clock::now();
        times_correct.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    for (std::size_t i = 0; i < N_REPLICA; i++) {
        if (qs_copy_impl[i] != qs_copy_correct[i]) {
            std::cerr << "Error in " << i << " th seq in " << name << ": reversed qual does not match original qual."
                      << std::endl;
            std::cerr << "ORIG: " << vec2str(qs[i]) << std::endl;
            std::cerr << "CORR: " << vec2str(qs_copy_correct[i]) << std::endl;
            std::cerr << "IMPL: " << vec2str(qs_copy_impl[i]) << std::endl;
            std::abort();
        }
    }
    std::cout << "IMPL: " << std::setw(25) << (name + ": ") << describe(times_impl) << std::endl;
    std::cout << "CORR: " << std::setw(25) << (name + ": ") << describe(times_correct) << std::endl;
}

template <typename T> void bench(const int rlen)
{
    std::size_t run_times = 0;

    run_times = BYTES_TO_PROCESS / rlen / sizeof(T);

    std::cout << ">run_times: " << run_times << " rlen: " << rlen << std::endl;

    std::vector<std::vector<T>> qs;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution dis(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());

    for (std::size_t i = 0; i < N_REPLICA; i++) {
        std::vector<T> q;
        q.reserve(rlen);
        for (int j = 0; j < rlen; j++) {
            q.emplace_back(dis(gen));
        }
        qs.emplace_back(q);
    }
    // Select the correct function ptr
    if constexpr (sizeof(T) == 8) {
        bech_impl(qReverse_8, "int-" + std::to_string(sizeof(T)), run_times, rlen, qs);
        return;
    }
    if constexpr (sizeof(T) == 4) {
        bech_impl(qReverse_4, "int-" + std::to_string(sizeof(T)), run_times, rlen, qs);
        return;
    }
    if constexpr (sizeof(T) == 2) {
        bech_impl(qReverse_2, "int-" + std::to_string(sizeof(T)), run_times, rlen, qs);
        return;
    }
    if constexpr (sizeof(T) == 1) {
        bech_impl(qReverse_1, "int-" + std::to_string(sizeof(T)), run_times, rlen, qs);
        return;
    }
}

} // namespace
int main()
{
    bench<int8_t>(36);
    bench<int8_t>(100);
    bench<int8_t>(150);
    bench<int8_t>(200);
    bench<int8_t>(300);
    bench<int8_t>(1000);
    bench<int8_t>(10000);

    bench<int16_t>(36);
    bench<int16_t>(100);
    bench<int16_t>(150);
    bench<int16_t>(200);
    bench<int16_t>(300);
    bench<int16_t>(1000);
    bench<int16_t>(10000);

    bench<int32_t>(36);
    bench<int32_t>(100);
    bench<int32_t>(150);
    bench<int32_t>(200);
    bench<int32_t>(300);
    bench<int32_t>(1000);
    bench<int32_t>(10000);

    bench<int64_t>(36);
    bench<int64_t>(100);
    bench<int64_t>(150);
    bench<int64_t>(200);
    bench<int64_t>(300);
    bench<int64_t>(1000);
    bench<int64_t>(10000);

    return 0;
}
