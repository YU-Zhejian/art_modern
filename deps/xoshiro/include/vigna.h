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

struct xoshiro_4x64_plus {
    uint64_t s[4];

    uint64_t next(void)
    {
        const uint64_t result = s[0] + s[3];

        const uint64_t t = s[1] << 17;

        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];

        s[2] ^= t;

        s[3] = rotl(s[3], 45);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^128 calls to next(); it can be used to generate 2^128
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
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
       2^192 calls to next(); it can be used to generate 2^64 starting points,
       from each of which jump() will generate 2^64 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241,
                                             0x39109bb02acbe635};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
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

struct xoshiro_4x64_plus_plus {
    uint64_t s[4];

    uint64_t next(void)
    {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^128 calls to next(); it can be used to generate 2^128
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
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
       2^192 calls to next(); it can be used to generate 2^64 starting points,
       from each of which jump() will generate 2^64 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241,
                                             0x39109bb02acbe635};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
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

struct xoshiro_4x64_star_star {
    uint64_t s[4];

    uint64_t next(void)
    {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;

        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^128 calls to next(); it can be used to generate 2^128
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
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
       2^192 calls to next(); it can be used to generate 2^64 starting points,
       from each of which jump() will generate 2^64 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241,
                                             0x39109bb02acbe635};

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
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

struct xoshiro_8x64_plus {
    uint64_t s[8];

    uint64_t next(void)
    {
        const uint64_t result = s[0] + s[2];

        const uint64_t t = s[1] << 11;
        s[2] ^= s[0];
        s[5] ^= s[1];
        s[1] ^= s[2];
        s[7] ^= s[3];
        s[3] ^= s[4];
        s[4] ^= s[5];
        s[0] ^= s[6];
        s[6] ^= s[7];
        s[6] ^= t;
        s[7] = rotl(s[7], 21);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^256 calls to next(); it can be used to generate 2^256
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x33ed89b6e7a353f9, 0x760083d7955323be, 0x2837f2fbb5f22fae, 0x4b8c5674d309511c,
                                        0xb11ac47a7ba28c25, 0xf1be7667092bcc1c, 0x53851efdb6df0aaf, 0x1ebbc8b23eaf25db};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b)
                    for (unsigned w = 0; w < sizeof s / sizeof *s; w++) t[w] ^= s[w];
                next();
            }

        memcpy(s, t, sizeof s);
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^384 calls to next(); it can be used to generate 2^128 starting points,
       from each of which jump() will generate 2^128 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x11467fef8f921d28, 0xa2a819f2e79c8ea8, 0xa8299fc284b3959a,
                                             0xb4d347340ca63ee1, 0x1cb0940bedbff6ce, 0xd956c5c4fa1f8e17,
                                             0x915e38fd4eda93bc, 0x5b3ccdfa5d7daca5};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b)
                    for (unsigned w = 0; w < sizeof s / sizeof *s; w++) t[w] ^= s[w];
                next();
            }

        memcpy(s, t, sizeof s);
    }
};

struct xoshiro_8x64_plus_plus {
    uint64_t s[8];

    uint64_t next(void)
    {
        const uint64_t result = rotl(s[0] + s[2], 17) + s[2];

        const uint64_t t = s[1] << 11;
        s[2] ^= s[0];
        s[5] ^= s[1];
        s[1] ^= s[2];
        s[7] ^= s[3];
        s[3] ^= s[4];
        s[4] ^= s[5];
        s[0] ^= s[6];
        s[6] ^= s[7];
        s[6] ^= t;
        s[7] = rotl(s[7], 21);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^256 calls to next(); it can be used to generate 2^256
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x33ed89b6e7a353f9, 0x760083d7955323be, 0x2837f2fbb5f22fae, 0x4b8c5674d309511c,
                                        0xb11ac47a7ba28c25, 0xf1be7667092bcc1c, 0x53851efdb6df0aaf, 0x1ebbc8b23eaf25db};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b)
                    for (unsigned w = 0; w < sizeof s / sizeof *s; w++) t[w] ^= s[w];
                next();
            }

        memcpy(s, t, sizeof s);
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^384 calls to next(); it can be used to generate 2^128 starting points,
       from each of which jump() will generate 2^128 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x11467fef8f921d28, 0xa2a819f2e79c8ea8, 0xa8299fc284b3959a,
                                             0xb4d347340ca63ee1, 0x1cb0940bedbff6ce, 0xd956c5c4fa1f8e17,
                                             0x915e38fd4eda93bc, 0x5b3ccdfa5d7daca5};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b)
                    for (unsigned w = 0; w < sizeof s / sizeof *s; w++) t[w] ^= s[w];
                next();
            }

        memcpy(s, t, sizeof s);
    }
};

struct xoshiro_8x64_star_star {
    uint64_t s[8];

    uint64_t next(void)
    {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;

        const uint64_t t = s[1] << 11;
        s[2] ^= s[0];
        s[5] ^= s[1];
        s[1] ^= s[2];
        s[7] ^= s[3];
        s[3] ^= s[4];
        s[4] ^= s[5];
        s[0] ^= s[6];
        s[6] ^= s[7];
        s[6] ^= t;
        s[7] = rotl(s[7], 21);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^256 calls to next(); it can be used to generate 2^256
       non-overlapping subsequences for parallel computations. */

    void jump(void)
    {
        static const uint64_t JUMP[] = {0x33ed89b6e7a353f9, 0x760083d7955323be, 0x2837f2fbb5f22fae, 0x4b8c5674d309511c,
                                        0xb11ac47a7ba28c25, 0xf1be7667092bcc1c, 0x53851efdb6df0aaf, 0x1ebbc8b23eaf25db};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b)
                    for (unsigned w = 0; w < sizeof s / sizeof *s; w++) t[w] ^= s[w];
                next();
            }

        memcpy(s, t, sizeof s);
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^384 calls to next(); it can be used to generate 2^128 starting points,
       from each of which jump() will generate 2^128 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {0x11467fef8f921d28, 0xa2a819f2e79c8ea8, 0xa8299fc284b3959a,
                                             0xb4d347340ca63ee1, 0x1cb0940bedbff6ce, 0xd956c5c4fa1f8e17,
                                             0x915e38fd4eda93bc, 0x5b3ccdfa5d7daca5};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b)
                    for (unsigned w = 0; w < sizeof s / sizeof *s; w++) t[w] ^= s[w];
                next();
            }

        memcpy(s, t, sizeof s);
    }
};

struct xoroshiro_16x64_plus_plus {
    unsigned p = 15; // Original had p initialized to 0 but that is an off-by-one thing
    uint64_t s[16];

    uint64_t next(void)
    {
        const unsigned q = p;
        const uint64_t s0 = s[p = (p + 1) & 15];
        uint64_t       s15 = s[q];
        const uint64_t result = rotl(s0 + s15, 23) + s15;

        s15 ^= s0;
        s[q] = rotl(s0, 25) ^ s15 ^ (s15 << 27);
        s[p] = rotl(s15, 36);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^512 calls to next(); it can be used to generate 2^512
       non-overlapping subsequences for parallel computations. */

    void jump()
    {
        static const uint64_t JUMP[] = {0x931197d8e3177f17, 0xb59422e0b9138c5f, 0xf06a6afb49d668bb, 0xacb8a6412c8a1401,
                                        0x12304ec85f0b3468, 0xb7dfe7079209891e, 0x405b7eec77d9eb14, 0x34ead68280c44e4a,
                                        0xe0e4ba3e0ac9e366, 0x8f46eda8348905b7, 0x328bf4dbad90d6ff, 0xc8fd6fb31c9effc3,
                                        0xe899d452d4b67652, 0x45f387286ade3205, 0x03864f454a8920bd, 0xa68fa28725b1b384};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b)
                    for (unsigned j = 0; j < sizeof s / sizeof *s; j++) t[j] ^= s[(j + p) & (sizeof s / sizeof *s - 1)];
                next();
            }

        for (unsigned i = 0; i < sizeof s / sizeof *s; i++) { s[(i + p) & (sizeof s / sizeof *s - 1)] = t[i]; }
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^768 calls to next(); it can be used to generate 2^256 starting points,
       from each of which jump() will generate 2^256 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {
            0x7374156360bbf00f, 0x4630c2efa3b3c1f6, 0x6654183a892786b1, 0x94f7bfcbfb0f1661,
            0x27d8243d3d13eb2d, 0x9701730f3dfb300f, 0x2f293baae6f604ad, 0xa661831cb60cd8b6,
            0x68280c77d9fe008c, 0x50554160f5ba9459, 0x2fc20b17ec7b2a9a, 0x49189bbdc8ec9f8f,
            0x92a65bca41852cc1, 0xf46820dd0509c12a, 0x52b00c35fbf92185, 0x1e5b3b7f589e03c1};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b)
                    for (unsigned j = 0; j < sizeof s / sizeof *s; j++) t[j] ^= s[(j + p) & (sizeof s / sizeof *s - 1)];
                next();
            }

        for (unsigned i = 0; i < sizeof s / sizeof *s; i++) { s[(i + p) & (sizeof s / sizeof *s - 1)] = t[i]; }
    }
};

struct xoroshiro_16x64_star {
    unsigned p = 15; // Original had p initialized to 0 but that is an off-by-one thing
    uint64_t s[16];

    uint64_t next(void)
    {
        const unsigned q = p;
        const uint64_t s0 = s[p = (p + 1) & 15];
        uint64_t       s15 = s[q];
        const uint64_t result = s0 * 0x9e3779b97f4a7c13;

        s15 ^= s0;
        s[q] = rotl(s0, 25) ^ s15 ^ (s15 << 27);
        s[p] = rotl(s15, 36);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^512 calls to next(); it can be used to generate 2^512
       non-overlapping subsequences for parallel computations. */

    void jump()
    {
        static const uint64_t JUMP[] = {0x931197d8e3177f17, 0xb59422e0b9138c5f, 0xf06a6afb49d668bb, 0xacb8a6412c8a1401,
                                        0x12304ec85f0b3468, 0xb7dfe7079209891e, 0x405b7eec77d9eb14, 0x34ead68280c44e4a,
                                        0xe0e4ba3e0ac9e366, 0x8f46eda8348905b7, 0x328bf4dbad90d6ff, 0xc8fd6fb31c9effc3,
                                        0xe899d452d4b67652, 0x45f387286ade3205, 0x03864f454a8920bd, 0xa68fa28725b1b384};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b)
                    for (unsigned j = 0; j < sizeof s / sizeof *s; j++) t[j] ^= s[(j + p) & (sizeof s / sizeof *s - 1)];
                next();
            }

        for (unsigned i = 0; i < sizeof s / sizeof *s; i++) { s[(i + p) & (sizeof s / sizeof *s - 1)] = t[i]; }
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^768 calls to next(); it can be used to generate 2^256 starting points,
       from each of which jump() will generate 2^256 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {
            0x7374156360bbf00f, 0x4630c2efa3b3c1f6, 0x6654183a892786b1, 0x94f7bfcbfb0f1661,
            0x27d8243d3d13eb2d, 0x9701730f3dfb300f, 0x2f293baae6f604ad, 0xa661831cb60cd8b6,
            0x68280c77d9fe008c, 0x50554160f5ba9459, 0x2fc20b17ec7b2a9a, 0x49189bbdc8ec9f8f,
            0x92a65bca41852cc1, 0xf46820dd0509c12a, 0x52b00c35fbf92185, 0x1e5b3b7f589e03c1};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b)
                    for (unsigned j = 0; j < sizeof s / sizeof *s; j++) t[j] ^= s[(j + p) & (sizeof s / sizeof *s - 1)];
                next();
            }

        for (unsigned i = 0; i < sizeof s / sizeof *s; i++) { s[(i + p) & (sizeof s / sizeof *s - 1)] = t[i]; }
    }
};

struct xoroshiro_16x64_star_star {
    unsigned p = 15; // Original had p initialized to 0 but that is an off-by-one thing
    uint64_t s[16];

    uint64_t next(void)
    {
        const unsigned q = p;
        const uint64_t s0 = s[p = (p + 1) & 15];
        uint64_t       s15 = s[q];
        const uint64_t result = rotl(s0 * 5, 7) * 9;

        s15 ^= s0;
        s[q] = rotl(s0, 25) ^ s15 ^ (s15 << 27);
        s[p] = rotl(s15, 36);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^512 calls to next(); it can be used to generate 2^512
       non-overlapping subsequences for parallel computations. */

    void jump()
    {
        static const uint64_t JUMP[] = {0x931197d8e3177f17, 0xb59422e0b9138c5f, 0xf06a6afb49d668bb, 0xacb8a6412c8a1401,
                                        0x12304ec85f0b3468, 0xb7dfe7079209891e, 0x405b7eec77d9eb14, 0x34ead68280c44e4a,
                                        0xe0e4ba3e0ac9e366, 0x8f46eda8348905b7, 0x328bf4dbad90d6ff, 0xc8fd6fb31c9effc3,
                                        0xe899d452d4b67652, 0x45f387286ade3205, 0x03864f454a8920bd, 0xa68fa28725b1b384};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b)
                    for (unsigned j = 0; j < sizeof s / sizeof *s; j++) t[j] ^= s[(j + p) & (sizeof s / sizeof *s - 1)];
                next();
            }

        for (unsigned i = 0; i < sizeof s / sizeof *s; i++) { s[(i + p) & (sizeof s / sizeof *s - 1)] = t[i]; }
    }

    /* This is the long-jump function for the generator. It is equivalent to
       2^768 calls to next(); it can be used to generate 2^256 starting points,
       from each of which jump() will generate 2^256 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void)
    {
        static const uint64_t LONG_JUMP[] = {
            0x7374156360bbf00f, 0x4630c2efa3b3c1f6, 0x6654183a892786b1, 0x94f7bfcbfb0f1661,
            0x27d8243d3d13eb2d, 0x9701730f3dfb300f, 0x2f293baae6f604ad, 0xa661831cb60cd8b6,
            0x68280c77d9fe008c, 0x50554160f5ba9459, 0x2fc20b17ec7b2a9a, 0x49189bbdc8ec9f8f,
            0x92a65bca41852cc1, 0xf46820dd0509c12a, 0x52b00c35fbf92185, 0x1e5b3b7f589e03c1};

        uint64_t t[sizeof s / sizeof *s];
        memset(t, 0, sizeof t);
        for (unsigned i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b)
                    for (unsigned j = 0; j < sizeof s / sizeof *s; j++) t[j] ^= s[(j + p) & (sizeof s / sizeof *s - 1)];
                next();
            }

        for (unsigned i = 0; i < sizeof s / sizeof *s; i++) { s[(i + p) & (sizeof s / sizeof *s - 1)] = t[i]; }
    }
};

} // namespace old
