#ifndef RANDSTR_H_INCLUDED
#define RANDSTR_H_INCLUDED

#include "pcg32c_minimal.h"

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline void generate_randstr(char* out, size_t length, pcg32_random_t* rng)
{
    static const char charset[] = "0123456789"
                                  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                  "abcdefghijklmnopqrstuvwxyz";
    const size_t charset_size = sizeof(charset) - 1;
    for (size_t i = 0; i < length; i++) {
        uint32_t index = pcg32_boundedrand_r(rng, (uint32_t)charset_size);
        out[i] = charset[index];
    }
    out[length] = '\0';
}

static inline void generate_randbuf(uint8_t* out, size_t length, pcg32_random_t* rng)
{
    for (size_t i = 0; i < length; i++) {
        out[i] = (uint8_t)pcg32_boundedrand_r(rng, 256);
    }
}

#ifdef __cplusplus
}
#endif

#endif // RANDSTR_H_INCLUDED
