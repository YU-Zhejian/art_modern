/**
 * Copyright 2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Dtypes.hh"
#include "libam_support/utils/seq_utils.hh"

#include <htslib/sam.h>

#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace labw::art_modern; // NOLINT

constexpr std::size_t N_REPLICA = 50UL;

const std::string benchmark_using_bam = "/dev/null"; // Replace with a real BAM file path for actual benchmarking

int main()
{
    // Read the BAM file and collect all CIGARs
    samFile* in = CExceptionsProxy::assert_not_null(sam_open(benchmark_using_bam.c_str(), "r"), "HTSLib");
    bam_hdr_t* header = CExceptionsProxy::assert_not_null(sam_hdr_read(in), "HTSLib");
    bam1_t* b = bam_init1();
    std::vector<std::vector<am_cigar_t>> cigars;
    while (sam_read1(in, header, b) >= 0) {
        cigars.emplace_back();
        cigars.back().resize(b->core.n_cigar);
        std::memcpy(cigars.back().data(), bam_get_cigar(b), b->core.n_cigar * sizeof(am_cigar_t));
    }
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    std::cout << "Read " << format_with_commas(cigars.size()) << " reads." << std::endl;

    // Benchmark

    std::vector<std::size_t> times;
    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto& cigars : cigars) {
            volatile auto s = cigar_arr_to_str_old(cigars.data(), cigars.size()); // NOLINT
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << std::setw(25) << "old: " << describe(times) << std::endl;
    times.clear();
    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto& cigars : cigars) {
            volatile auto s = cigar_arr_to_str_optim(cigars.data(), cigars.size()); // NOLINT
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::cout << std::setw(25) << "new: " << describe(times) << std::endl;
}
