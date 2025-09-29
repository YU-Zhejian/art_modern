/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

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
    GslRngWrapper(const gsl_rng_type* t, const result_type seed)
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
