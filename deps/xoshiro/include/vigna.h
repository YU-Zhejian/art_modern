/// @file vigna.h
/// @brief Here are the xoshiro/xoroshiro "C" implementations from the author's website.
///
/// Each version is simply the original pasted into a struct with an appropriate name where the struct owns the state
/// array for the particular algorithm. Didn't really want to alter the original code in any material way--to get rid of
/// some signed/unsigned compiler comparison warnings we did change a few variables from @c int to @c unsigned. For the
/// same reason we added some parentheses suggested by a later version of GCC.
///
/// The only other change we made was in the Xoroshiro 1024 bit versions -- the originals default initialized the cycle
/// variable `p` to zero. To make those algorithms consistent with our Xoroshiro implementations we instead initialize
/// that same variable p  to 15. This means that on the first call we mix s[0] and s[15] and then on in a cycle from
/// there. The original starts by mixing s[1] and s[0] and on in the cycle from there -- a trivial off-by-one change.
///
/// @copyright The originals all carry the following copyright notice ...
/// Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)
/// To the extent possible under law, the author has dedicated all copyright and related and neighboring rights to this
/// software to the public domain worldwide. This software is distributed without any warranty.
/// See <http://creativecommons.org/publicdomain/zero/1.0/>.
#pragma once

#include <iostream>
#include <stdint.h>
#include <string.h>

namespace old {

static inline uint32_t
rotl(const uint32_t x, int k)
{
    return (x << k) | (x >> (32 - k));
}

static inline uint64_t
rotl(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

struct xoroshiro_2x32_star {
    uint32_t s[2];

    uint32_t next(void)
    {
        const uint32_t s0 = s[0];
        uint32_t       s1 = s[1];
        const uint32_t result = s0 * 0x9E3779BB;

        s1 ^= s0;
        s[0] = rotl(s0, 26) ^ s1 ^ (s1 << 9); // a, b
        s[1] = rotl(s1, 13);                  // c

        return result;
    }
};

struct xoroshiro_2x32_star_star {
    uint32_t s[2];
    uint32_t next(void)
    {
        const uint32_t s0 = s[0];
        uint32_t       s1 = s[1];
        const uint32_t result = rotl(s0 * 0x9E3779BB, 5) * 5;

        s1 ^= s0;
        s[0] = rotl(s0, 26) ^ s1 ^ (s1 << 9); // a, b
        s[1] = rotl(s1, 13);                  // c

        return result;
    }
};

struct xoshiro_4x32_plus {
    uint32_t s[4];

    uint32_t next(void)
    {
        const uint32_t result = s[0] + s[3];

        const uint32_t t = s[1] << 9;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 11);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^64 calls to next(); it can be used to generate 2^64
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint32_t JUMP[] = {0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b};

        uint32_t s0 = 0;
        uint32_t s1 = 0;
        uint32_t s2 = 0;
        uint32_t s3 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 32; b++) {
                if (JUMP[i] & UINT32_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^96 calls to next(); it can be used to generate 2^32 starting points,
       from each of which jump() will generate 2^32 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint32_t LONG_JUMP[] = {0xb523952e, 0x0b6f099f, 0xccf5a0ef, 0x1c580662};

        uint32_t s0 = 0;
        uint32_t s1 = 0;
        uint32_t s2 = 0;
        uint32_t s3 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 32; b++) {
                if (LONG_JUMP[i] & UINT32_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }
};

struct xoshiro_4x32_plus_plus {
    uint32_t s[4];

    uint32_t next(void)
    {
        const uint32_t result = rotl(s[0] + s[3], 7) + s[0];

        const uint32_t t = s[1] << 9;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 11);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^64 calls to next(); it can be used to generate 2^64
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint32_t JUMP[] = {0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b};

        uint32_t s0 = 0;
        uint32_t s1 = 0;
        uint32_t s2 = 0;
        uint32_t s3 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 32; b++) {
                if (JUMP[i] & UINT32_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^96 calls to next(); it can be used to generate 2^32 starting points,
       from each of which jump() will generate 2^32 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint32_t LONG_JUMP[] = {0xb523952e, 0x0b6f099f, 0xccf5a0ef, 0x1c580662};

        uint32_t s0 = 0;
        uint32_t s1 = 0;
        uint32_t s2 = 0;
        uint32_t s3 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 32; b++) {
                if (LONG_JUMP[i] & UINT32_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }
};

struct xoshiro_4x32_star_star {
    uint32_t s[4];

    uint32_t next(void)
    {
        const uint32_t result = rotl(s[1] * 5, 7) * 9;

        const uint32_t t = s[1] << 9;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 11);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^64 calls to next(); it can be used to generate 2^64
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint32_t JUMP[] = {0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b};

        uint32_t s0 = 0;
        uint32_t s1 = 0;
        uint32_t s2 = 0;
        uint32_t s3 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 32; b++) {
                if (JUMP[i] & UINT32_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^96 calls to next(); it can be used to generate 2^32 starting points,
       from each of which jump() will generate 2^32 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint32_t LONG_JUMP[] = {0xb523952e, 0x0b6f099f, 0xccf5a0ef, 0x1c580662};

        uint32_t s0 = 0;
        uint32_t s1 = 0;
        uint32_t s2 = 0;
        uint32_t s3 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 32; b++) {
                if (LONG_JUMP[i] & UINT32_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }
};

struct xoroshiro_2x64_plus {
    uint64_t s[2];

    uint64_t next(void)
    {
        const uint64_t s0 = s[0];
        uint64_t       s1 = s[1];
        const uint64_t result = s0 + s1;

        s1 ^= s0;
        s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        s[1] = rotl(s1, 37);                   // c

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^64 calls to next(); it can be used to generate 2^64
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0xdf900294d8f554a5, 0x170865df4b3201fc};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^96 calls to next(); it can be used to generate 2^32 starting points,
       from each of which jump() will generate 2^32 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
    }
};

struct xoroshiro_2x64_plus_plus {
    uint64_t s[2];

    uint64_t next(void)
    {
        const uint64_t s0 = s[0];
        uint64_t       s1 = s[1];
        const uint64_t result = rotl(s0 + s1, 17) + s0;

        s1 ^= s0;
        s[0] = rotl(s0, 49) ^ s1 ^ (s1 << 21); // a, b
        s[1] = rotl(s1, 28);                   // c

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^64 calls to next(); it can be used to generate 2^64
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x2bd7a6a6e99c2ddc, 0x0992ccaf6a6fca05};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^96 calls to next(); it can be used to generate 2^32 starting points,
       from each of which jump() will generate 2^32 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x360fd5f2cf8d5d99, 0x9c6e6877736c46e3};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
    }
};

struct xoroshiro_2x64_star_star {
    uint64_t s[2];

    uint64_t next(void)
    {
        const uint64_t s0 = s[0];
        uint64_t       s1 = s[1];
        const uint64_t result = rotl(s0 * 5, 7) * 9;

        s1 ^= s0;
        s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        s[1] = rotl(s1, 37);                   // c

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^64 calls to next(); it can be used to generate 2^64
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0xdf900294d8f554a5, 0x170865df4b3201fc};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^96 calls to next(); it can be used to generate 2^32 starting points,
       from each of which jump() will generate 2^32 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
    }
};

} // namespace old
