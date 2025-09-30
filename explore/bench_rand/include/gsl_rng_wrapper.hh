#pragma once

#include "class_utils.hh"

#include <absl/base/attributes.h>

#include <gsl/gsl_rng.h>

#include <string>

class GslRngWrapper {
public:
    DELETE_COPY_MOVE(GslRngWrapper)
    using result_type = decltype(gsl_rng_get(nullptr));
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
    ABSL_ATTRIBUTE_ALWAYS_INLINE result_type operator()() { return gsl_rng_get(r); }
    ABSL_ATTRIBUTE_ALWAYS_INLINE result_type min() { return gsl_rng_min(r); }
    ABSL_ATTRIBUTE_ALWAYS_INLINE result_type max() { return gsl_rng_max(r); }
    [[nodiscard]] std::string name() const { return gsl_rng_name(r); }

private:
    gsl_rng* r;

}; // namespace labw::art_modern
