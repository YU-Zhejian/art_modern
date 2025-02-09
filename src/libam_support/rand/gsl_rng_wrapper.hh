#pragma once

#include "libam_support/utils/class_macros_utils.hh"

#include <gsl/gsl_rng.h>

#include <string>

namespace labw::art_modern {
class GslRngWrapper {
public:
    DELETE_COPY(GslRngWrapper)
    DELETE_MOVE(GslRngWrapper)
    using result_type = unsigned long;
    explicit GslRngWrapper(const gsl_rng_type* t)
        : r(gsl_rng_alloc(t))
    {
    }
    GslRngWrapper(const gsl_rng_type* t, result_type seed)
        : GslRngWrapper(t)
    {
        gsl_rng_set(r, seed);
    }
    ~GslRngWrapper() { gsl_rng_free(r); }
    result_type operator()() { return gsl_rng_get(r); }
    result_type min() { return gsl_rng_min(r); }
    result_type max() { return gsl_rng_max(r); }
    [[nodiscard]] std::string name() const { return gsl_rng_name(r); }

private:
    gsl_rng* r;
};

} // namespace labw::art_modern
