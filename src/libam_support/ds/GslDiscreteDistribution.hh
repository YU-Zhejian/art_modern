/**
 * GSL Discrete Distribution
 *
 * Re-implemented in C++ from <randist/discrete.c> in GSL 2.8.
 * Licensed under GPL.
 */
#pragma once

#include <cmath>
#include <cstdlib>
#include <numeric>
#include <stack>
#include <vector>

namespace labw::art_modern {
template <typename FloatType> class GslDiscreteDistribution {
    size_t K;
    std::vector<FloatType> F;
    std::vector<size_t> A;

public:
    GslDiscreteDistribution()
        : K(0) {};
    explicit GslDiscreteDistribution(const std::vector<FloatType>& prob_array)
        : K(prob_array.size())
        , F(std::vector<FloatType>(prob_array.size()))
        , A(std::vector<size_t>(prob_array.size()))
    {
        std::stack<size_t> Bigs;
        std::stack<size_t> Smalls;
        FloatType const pTotal = std::accumulate(prob_array.begin(), prob_array.end(), 0.0);
        FloatType mean = 0;
        std::vector<FloatType> E(K);

        for (size_t k = 0; k < K; ++k) {
            E[k] = prob_array[k] / pTotal;
        }

        /* Now create the Bigs and the Smalls */
        mean = 1.0 / K;

        for (size_t k = 0; k < K; ++k) {
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
        for (size_t k = 0; k < K; ++k) {
            F[k] += k;
            F[k] /= K;
        }
#endif
    }

    /**
     * @param u An uniform random number in [0,1)
     * @return
     */
    size_t operator()(FloatType u) const
    {
        size_t c = 0;
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
};
} // namespace labw::art_modern
