#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <chrono>
#include <cstdint>

/**
 * Minimalist PCG32 implementation in C.
 *
 * See: <https://github.com/imneme/pcg-c-basic>
 *
 * *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
 * Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
 */
class pcg32_c {
public:
    pcg32_c(uint64_t state, uint64_t inc)
        : state_(state)
        , inc_(inc) {};

    pcg32_c()
        : pcg32_c { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL } {};
    DEFAULT_DESTRUCTOR(pcg32_c)
    DELETE_COPY(pcg32_c)
    DELETE_MOVE(pcg32_c)
    using result_type = uint32_t;
    constexpr result_type min() const noexcept { return std::numeric_limits<result_type>::min(); }
    constexpr result_type max() const noexcept { return std::numeric_limits<result_type>::max(); }

    uint32_t operator()()
    {
        uint64_t const oldstate = state_;
        // Advance internal state
        state_ = oldstate * 6364136223846793005ULL + (inc_ | 1);
        // Calculate output function (XSH RR), uses old state for max ILP
        uint32_t const xorshifted = ((oldstate >> 18U) ^ oldstate) >> 27U;
        uint32_t const rot = oldstate >> 59U;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

    void seed(const uint64_t initstate, const uint64_t initseq)
    {
        state_ = 0U;
        inc_ = (initseq << 1U) | 1U;
        (void)(*this)();
        state_ += initstate;
        (void)(*this)();
    }

    void seed()
    {
        auto now = std::chrono::system_clock::now().time_since_epoch().count();
        auto addr = reinterpret_cast<intptr_t>(this);
        seed(now, addr);
    }

private:
    uint64_t state_;
    uint64_t inc_;
};
