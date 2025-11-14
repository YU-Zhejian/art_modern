#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <cstdint>
#include <limits>

namespace labw::art_modern {

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
    using result_type = uint32_t;
    pcg32_c(uint64_t state, uint64_t inc);
    explicit pcg32_c(uint64_t seed_value);
    pcg32_c();

    DEFAULT_DESTRUCTOR(pcg32_c)
    DELETE_COPY(pcg32_c)
    DELETE_MOVE(pcg32_c)

    static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

    result_type operator()();

    void seed(uint64_t initstate, uint64_t initseq);
    void seed(uint64_t seed_value);

private:
    uint64_t state_;
    uint64_t inc_;
};
} // namespace labw::art_modern
