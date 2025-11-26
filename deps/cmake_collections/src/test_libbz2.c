#include "randstr.h"

#include <bzlib.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define THIS_BZ_VERBOSITY 3

int perform(uint8_t* randbuf, size_t randbuf_len)
{
    int retv;
    unsigned int max_compressed_size = randbuf_len + (randbuf_len / 100) + 600;
    uint8_t* compressed_buf = (uint8_t*)malloc(max_compressed_size);
    unsigned int compressed_size = max_compressed_size;
    if (compressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }
    retv = BZ2_bzBuffToBuffCompress(
        (char*)compressed_buf, &compressed_size, (char*)randbuf, (unsigned int)randbuf_len, 1, THIS_BZ_VERBOSITY, 0);
    if (retv != BZ_OK) {
        fprintf(stderr, "Compression failed with error code %d\n", retv);
        free(compressed_buf);
        return 1;
    }
    unsigned int decompressed_size = (unsigned int)randbuf_len;
    uint8_t* decompressed_buf = (uint8_t*)malloc(decompressed_size);
    if (decompressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(compressed_buf);
        return 1;
    }
    retv = BZ2_bzBuffToBuffDecompress(
        (char*)decompressed_buf, &decompressed_size, (char*)compressed_buf, compressed_size, 0, THIS_BZ_VERBOSITY);
    if (retv != BZ_OK) {
        fprintf(stderr, "Decompression failed with error code %d\n", retv);
        free(compressed_buf);
        free(decompressed_buf);
        return 1;
    }

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
    printf("libbz2 version: %s\n", BZ2_bzlibVersion());
    const size_t randbuf_len = 4ULL * 1024;
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
