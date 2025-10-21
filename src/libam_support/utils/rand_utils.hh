#pragma once

#include <cstdint>

namespace labw::art_modern {
/**
 * Generate a random seed using current time, thread id, MPI rank and random device.
 * @return The generated random seed.
 */
std::uint64_t rand_seed();
} // namespace labw::art_modern
