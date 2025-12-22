#pragma once
#include "art_modern_config.h" // NOLINT: For WITH_MPI
#include "libam_support/Dtypes.h"

#include "libam_support/utils/class_macros_utils.hh"

#include <cstdlib>
#include <vector>
namespace labw::art_modern {
class SeedAlloc {
public:
    SeedAlloc() = default;
    ~SeedAlloc() = default;
    DELETE_COPY(SeedAlloc)
    DELETE_MOVE(SeedAlloc)

    static constexpr int MASTER_SEED_MPI_TAG = 1234;

    void run_seedalloc(bool ignore_seed_from_cmdline, am_rand_seed_t seed_from_cmdline, std::size_t n_threads);

    [[nodiscard]] am_rand_seed_t seed(std::size_t thread_idx) const;

private:
    am_rand_seed_t master_seed_ { 0 };
    am_rand_seed_t this_process_master_seed_ { 0 };
    std::vector<am_rand_seed_t> allocated_seeds_;
};
} // namespace  labw::art_modern
