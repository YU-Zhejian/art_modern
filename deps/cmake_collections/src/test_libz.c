#include "randstr.h"

#include <zlib.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int perform(uint8_t* randbuf, size_t randbuf_len)
{
    int retv;
    // Estimate the maximum compressed size
    size_t max_compressed_size = compressBound(randbuf_len);
    uint8_t* compressed_buf = (uint8_t*)malloc(max_compressed_size);
    uLongf compressed_size = (uLongf)max_compressed_size;
    if (compressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    retv = compress((Bytef*)compressed_buf, &compressed_size, (Bytef*)randbuf, (uLongf)randbuf_len);
    if (retv != Z_OK) {
        fprintf(stderr, "Compression failed with error code %d\n", retv);
        free(compressed_buf);
        return 1;
    }
    uLongf decompressed_size = (uLongf)randbuf_len;
    uint8_t* decompressed_buf = (uint8_t*)malloc(decompressed_size);
    if (decompressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(compressed_buf);
        return 1;
    }
    retv = uncompress(decompressed_buf, &decompressed_size, compressed_buf, compressed_size);
    if (retv != Z_OK) {
        fprintf(stderr, "Decompression failed with error code %d\n", retv);
        free(compressed_buf);
        free(decompressed_buf);
        return 1;
    }
    // Verify that the decompressed data matches the original data
    if (decompressed_size != randbuf_len || memcmp(randbuf, decompressed_buf, randbuf_len) != 0) {
        fprintf(stderr, "Decompressed data does not match original data\n");
        free(compressed_buf);
        free(decompressed_buf);
        return 1;
    }
    free(compressed_buf);
    free(decompressed_buf);
    return 0;
}

int main(void)
{
    printf("zlib version: %s\n", ZLIB_VERSION);
    const size_t randbuf_len = 4ULL * 1024; // 4kB
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (uint64_t)time(NULL), (uint64_t)(uintptr_t)&rng);
    uint8_t* randbuf = (uint8_t*)malloc(randbuf_len);
    if (randbuf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }
    generate_randbuf(randbuf, randbuf_len, &rng);
    if (perform(randbuf, randbuf_len) != 0) {
        free(randbuf);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
