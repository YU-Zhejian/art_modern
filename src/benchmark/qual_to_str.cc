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
#include "libam_support/seq/qual_to_str.h"
#include "libam_support/seq/str_to_qual.h"
#include "libam_support/utils/seq_utils.hh"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace labw::art_modern; // NOLINT

constexpr std::size_t N_REPLICA = 2000UL;
constexpr std::size_t  BASES_TO_PROCESS = 20ULL * M_SIZE;

namespace {
    [[maybe_unused]] void qual_to_str_foreach(const am_qual_t* qual,char* str,  const size_t qlen)
    {
        std::memcpy(str, qual, qlen);
        std::for_each(str, str + qlen, [](char& c) { c += PHRED_OFFSET; });
    }

    [[maybe_unused]] void str_to_qual_foreach(am_qual_t* qual, const char* str, const size_t qlen)
    {
        std::memcpy(qual, str, qlen);
        std::for_each(qual, qual + qlen, [](am_qual_t& c) { c -= PHRED_OFFSET; });
    }


    void bech_impl(
            void (*f_qual_to_str)(const signed char*, char* str, std::size_t),
    void (*f_str_to_qual)(am_qual_t* qual, const char* str, size_t qlen),
    const std::string& name,
    const std::size_t run_times,
    const int rlen,
    const std::vector<std::vector<am_qual_t>>& qs
    )
{
    std::vector<am_qual_t> q_regen;
    q_regen.resize(rlen);
    std::vector<std::size_t> times;
    std::string s;
    s.resize(rlen);
//    std::string correct_s;
//    correct_s.resize(rlen);
    // f_qual_to_str(q.data(), correct_s.data(), rlen);
    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < run_times; i++) {
            f_qual_to_str(qs[j].data(), s.data(), rlen);
            f_str_to_qual(q_regen.data(), s.c_str(), rlen);
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
//    if (s != correct_s) {
//        std::cerr << "Error in " << name << ": regenerated string does not match original string." << std::endl;
//        std::abort();
//    }
//    if (q_regen != q) {
//        std::cerr << "Error in " << name << ": regenerated qual does not match original qual." << std::endl;
//        std::abort();
//    }
    std::cout << std::setw(25) << (name + ": ") << describe(times) << std::endl;
}

void bench(const int rlen)
{
    std::size_t const run_times = BASES_TO_PROCESS / rlen;
    std::cout << ">run_times: " << run_times << " rlen: " << rlen << std::endl;
    std::vector< std::vector<am_qual_t>> qs;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution dis(0, 40);

    for (std::size_t i = 0; i < N_REPLICA; i++) {
        std::vector<am_qual_t> q;
        q.reserve(rlen);
        for (int j = 0; j < rlen;j++) {
            q.emplace_back(dis(gen));
        }
        qs.emplace_back(q);
    }

    bech_impl(qual_to_str_for_loop, str_to_qual_for_loop, "Scala (for loop)", run_times, rlen, qs);
#ifdef __MMX__
    bech_impl(qual_to_str_mmx, str_to_qual_mmx, "MMX", run_times, rlen, qs);
#endif
#ifdef  __SSE2__
    bech_impl(qual_to_str_sse2, str_to_qual_sse2, "SSE2", run_times, rlen, qs);
#endif
#ifdef  __AVX2__
    bech_impl(qual_to_str_avx2, str_to_qual_avx2, "AVX2", run_times, rlen, qs);
#endif
}

} // namespace
int main()
{
    bench(36);
    bench(100);
    bench(150);
    bench(200);
    bench(300);
    bench(1000);
    bench(10000);
    return 0;
}
