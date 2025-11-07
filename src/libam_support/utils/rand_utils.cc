#include "libam_support/utils/rand_utils.hh"

#include "libam_support/utils/hash_utils.h"
#include "libam_support/utils/mpi_utils.hh"

#include <chrono>
#include <cstdint>
#include <random>
#include <thread>

namespace labw::art_modern {
std::uint64_t rand_seed()
{
    const std::uint64_t now
        = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
              .count();
    const std::uint64_t thread_id = std::hash<std::thread::id>()(std::this_thread::get_id());
    const std::uint64_t rank = have_mpi() ? mpi_rank() : 0;
    std::random_device rd;
    const std::uint64_t entropy = rd();
    std::uint64_t seed = 0;
    am_boost_hash_combine_impl64(seed, now);
    am_boost_hash_combine_impl64(seed, thread_id);
    am_boost_hash_combine_impl64(seed, rank);
    am_boost_hash_combine_impl64(seed, entropy);
    return seed;
}
} // namespace labw::art_modern
