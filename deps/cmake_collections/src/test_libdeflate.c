#include "randstr.h"

#include <libdeflate.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int perform(uint8_t* randbuf, size_t randbuf_len)
{
    int retv;
    struct libdeflate_compressor* c = libdeflate_alloc_compressor(7);
    struct libdeflate_decompressor* d = libdeflate_alloc_decompressor();
    if (c == NULL || d == NULL) {
        fprintf(stderr, "Failed to allocate compressor or decompressor\n");
        return 1;
    }

    // Estimate the maximum compressed size
    size_t max_compressed_size = libdeflate_gzip_compress_bound(c, randbuf_len);
    uint8_t* compressed_buf = (uint8_t*)malloc(max_compressed_size);
    if (compressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }
    size_t compressed_size = libdeflate_gzip_compress(c, randbuf, randbuf_len, compressed_buf, max_compressed_size);
    if (compressed_size == 0) {
        fprintf(stderr, "Compression failed\n");
        free(compressed_buf);
        return 1;
    }
    size_t decompressed_size = randbuf_len;
    uint8_t* decompressed_buf = (uint8_t*)malloc(decompressed_size);
    if (decompressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(compressed_buf);
        return 1;
    }
    retv = libdeflate_gzip_decompress(d, compressed_buf, compressed_size, decompressed_buf, decompressed_size, &decompressed_size);
    if (retv != LIBDEFLATE_SUCCESS) {
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
    // Why we do not need to free c and d?
    return 0;
}

int main(void)
{
    printf("libdeflate version: %s\n", LIBDEFLATE_VERSION_STRING);
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
