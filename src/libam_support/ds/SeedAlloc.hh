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
#pragma once
#include "art_modern_config.h" // NOLINT: For WITH_MPI
#include "libam_support/Dtypes.h"

#include "libam_support/ds/pcg_32_c.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstdlib>
#include <memory>
#include <vector>

namespace labw::art_modern {
class SeedAlloc {
public:
    SeedAlloc() = default;
    ~SeedAlloc() = default;
    DELETE_COPY(SeedAlloc)
    DELETE_MOVE(SeedAlloc)

    static constexpr int MASTER_SEED_MPI_TAG = 1234;

    void run_seedalloc(am_rand_seed_t seed);

    [[nodiscard]] am_rand_seed_t nextseed() const;

private:
    am_rand_seed_t master_seed_ { 0 };
    am_rand_seed_t this_process_master_seed_ { 0 };
    std::vector<am_rand_seed_t> allocated_seeds_;

    std::unique_ptr<pcg32_c> pcg_rng_;
};
} // namespace  labw::art_modern
