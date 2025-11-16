/**
 * GSL Discrete Distribution
 *
 * Re-implemented in C++ from <randist/discrete.c> in GSL 2.8.
 * Licensed under GPL.
 *
 * Original license as follows:
 *
 * 1996, 1997, 1998, 1999, 2000, 2007, 2009 James Theiler, Brian Gough
 * 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#pragma once

#include <cstdlib>
#include <numeric> // for std::accumulate
#include <stack>
#include <vector>

namespace labw::art_modern {
template <typename FloatType> class GslDiscreteDistribution {
    std::size_t K;
    std::vector<FloatType> F;
    std::vector<std::size_t> A;

public:
    GslDiscreteDistribution();
    explicit GslDiscreteDistribution(const std::vector<FloatType>& prob_array);

    /**
     * @param u A uniform random number in [0,1)
     * @return
     */
    std::size_t operator()(FloatType u) const;
};

template <typename FloatType>
GslDiscreteDistribution<FloatType>::GslDiscreteDistribution()
    : K(0)
{
}

template <typename FloatType>
GslDiscreteDistribution<FloatType>::GslDiscreteDistribution(const std::vector<FloatType>& prob_array)
    : K(prob_array.size())
    , F(std::vector<FloatType>(prob_array.size()))
    , A(std::vector<std::size_t>(prob_array.size()))
{
    std::stack<std::size_t> Bigs;
    std::stack<std::size_t> Smalls;
    FloatType const pTotal = std::accumulate(prob_array.begin(), prob_array.end(), 0.0);
    FloatType mean = 1.0 / K;
    std::vector<FloatType> E(K);

    for (std::size_t k = 0; k < K; ++k) {
        E[k] = prob_array[k] / pTotal;
    }

    /* Now create the Bigs and the Smalls */
    for (std::size_t k = 0; k < K; ++k) {
        auto& Dest = E[k] > mean ? Bigs : Smalls;
        Dest.push(k);
    }

    /* Now work through the smalls */
    while (!Smalls.empty()) {
        auto s = Smalls.top();
        Smalls.pop();
        if (Bigs.empty()) {
            A[s] = s;
            F[s] = 1.0;
            continue;
        }
        auto b = Bigs.top();
        Bigs.pop();
        A[s] = b;
        F[s] = K * E[s];
        FloatType const d = mean - E[s];
        E[s] += d; /* now E[s] == mean */
        E[b] -= d;
        if (E[b] < mean) {
            Smalls.push(b); /* no longer big, join ranks of the small */
        } else if (E[b] > mean) {
            Bigs.push(b); /* still big, put it back where you found it */
        } else {
            /* E[b]==mean implies it is finished too */
            A[b] = b;
            F[b] = 1.0;
        }
    }
    while (!Bigs.empty()) {
        auto b = Bigs.top();
        Bigs.pop();
        A[b] = b;
        F[b] = 1.0;
    }
#if (1) // KNUTH_CONVENTION
    /* For convenience, set F'[k]=(k+F[k])/K */
    /* This saves some arithmetic in gsl_ran_discrete(); I find that
     * it doesn't actually make much difference.
     */
    for (std::size_t k = 0; k < K; ++k) {
        F[k] += k;
        F[k] /= K;
    }
#endif
}

template <typename FloatType> std::size_t GslDiscreteDistribution<FloatType>::operator()(FloatType u) const
{
    std::size_t c = 0;
#if (1) // KNUTH_CONVENTION
    c = u * K;
#else
    u *= K;
    c = u;
    u -= c;
#endif
    FloatType const f = F[c];
    if (f == 1.0 || u < f) {
        return c;
    }
    return A[c];
}

} // namespace labw::art_modern
