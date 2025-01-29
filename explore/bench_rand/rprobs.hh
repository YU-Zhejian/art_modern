#pragma once

#include <cstddef>
#include <cstdint>
constexpr std::size_t N_BASES = (1 << 10);
constexpr std::size_t N_TIMES = 20 * (1 << 10);
constexpr std::size_t a = 0;
constexpr std::size_t b = 1000;
static uint64_t seed()
{
    return 0;
    //    return
    //    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
    //               .count()
    //        * static_cast<uint64_t>(std::hash<std::thread::id>()(std::this_thread::get_id()));
}