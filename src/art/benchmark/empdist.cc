/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
#include "art/lib/Empdist.hh"
#include "art/lib/ArtConstants.hh"
#include "art/lib/Rprob.hh"

#include "benchmark/benchmark_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/rand_utils.hh"

#include <boost/log/trivial.hpp>

#include <chrono>
#include <cstdlib>
#include <vector>

using namespace labw::art_modern;
constexpr std::size_t N_REPLICA = 500UL;
constexpr std::size_t N_TIMES = M_SIZE * static_cast<std::size_t>(10);

int main()
{
    Empdist const empdist { DEFAULT_ERR_PROFILE, false, true, false };

    Rprob rp(200.0, 10.0, 10, 150, rand_seed());
    std::vector<am_qual_t> qual(150);
    std::vector<std::size_t> times;

    for (std::size_t i = 0; i < N_REPLICA; i++) {
        const auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t j = 0; j < N_TIMES; j++) {
            empdist.get_read_qual(qual, rp, true);
            empdist.get_read_qual(qual, rp, false);
        }
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = end - start;
        times.emplace_back(std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count());
    }
    BOOST_LOG_TRIVIAL(info) << "Empdist::get_read_qual  : " << describe(times, true);
    // Intel: Empdist::get_read_qual  : gmean: 35,650,559; mean/sd: 35,743,146/2,945,319
    // GCC:  Empdist::get_read_qual  : gmean: 57,546,368; mean/sd: 59,443,214/20,806,297
    return EXIT_SUCCESS;
}
