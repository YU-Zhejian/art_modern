//
// Created by yuzj on 11/2/25.
//

#include "libam_support/ds/pcg_32_c.hh"

#include <cstdint>
#include <limits>

namespace labw::art_modern {
pcg32_c::pcg32_c(uint64_t state, uint64_t inc)
    : state_(state)
    , inc_(inc) {};
pcg32_c::pcg32_c()
    : pcg32_c { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL } {};
constexpr pcg32_c::result_type pcg32_c::min() noexcept { return std::numeric_limits<result_type>::min(); }
constexpr pcg32_c::result_type pcg32_c::max() noexcept { return std::numeric_limits<result_type>::max(); }
pcg32_c::result_type pcg32_c::operator()()
{
    uint64_t const oldstate = state_;
    // Advance internal state
    state_ = oldstate * 6364136223846793005ULL + (inc_ | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t const xorshifted = ((oldstate >> 18U) ^ oldstate) >> 27U;
    uint32_t const rot = oldstate >> 59U;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
pcg32_c::pcg32_c(uint64_t seed_value)
    : pcg32_c()
{
    seed(seed_value);
}
void pcg32_c::seed(const uint64_t initstate, const uint64_t initseq)
{
    state_ = 0U;
    inc_ = (initseq << 1U) | 1U;
    (void)(*this)();
    state_ += initstate;
    (void)(*this)();
}

void pcg32_c::seed(const uint64_t seed_value) { seed(seed_value, reinterpret_cast<intptr_t>(this)); }

} // namespace labw::art_modern
