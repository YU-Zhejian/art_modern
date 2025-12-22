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
