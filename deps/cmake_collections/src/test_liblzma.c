#include "randstr.h"

#include <lzma.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int perform(uint8_t* randbuf, size_t randbuf_len)
{
    lzma_ret ret;

    /* --- Compress with XZ (XZ container/wrapper) --- */
    lzma_stream enc = LZMA_STREAM_INIT;
    ret = lzma_easy_encoder(&enc, 6 /* preset */, LZMA_CHECK_CRC64);
    if (ret != LZMA_OK) {
        fprintf(stderr, "lzma_easy_encoder failed: %d\n", ret);
        return 1;
    }

    size_t out_bound = lzma_stream_buffer_bound(randbuf_len);
    uint8_t* compressed_buf = (uint8_t*)malloc(out_bound);
    if (compressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        lzma_end(&enc);
        return 1;
    }

    enc.next_in = randbuf;
    enc.avail_in = randbuf_len;
    enc.next_out = compressed_buf;
    enc.avail_out = out_bound;

    /* run until stream end */
    for (;;) {
        ret = lzma_code(&enc, LZMA_FINISH);
        if (ret == LZMA_STREAM_END)
            break;
        if (ret != LZMA_OK) {
            fprintf(stderr, "Compression failed: %d\n", ret);
            free(compressed_buf);
            lzma_end(&enc);
            return 1;
        }
        /* If we loop here because avail_out became zero, that's an unexpected
           condition because we sized the buffer with lzma_stream_buffer_bound(). */
    }

    size_t compressed_size = (size_t)(enc.next_out - compressed_buf);
    lzma_end(&enc);

    /* --- Decompress XZ stream --- */
    lzma_stream dec = LZMA_STREAM_INIT;
    ret = lzma_stream_decoder(&dec, UINT64_MAX, LZMA_CONCATENATED | LZMA_FAIL_FAST);
    if (ret != LZMA_OK) {
        fprintf(stderr, "lzma_stream_decoder failed: %d\n", ret);
        free(compressed_buf);
        return 1;
    }

    uint8_t* decompressed_buf = (uint8_t*)malloc(randbuf_len);
    if (decompressed_buf == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(compressed_buf);
        lzma_end(&dec);
        return 1;
    }

    dec.next_in = compressed_buf;
    dec.avail_in = compressed_size;
    dec.next_out = decompressed_buf;
    dec.avail_out = randbuf_len;

    for (;;) {
        ret = lzma_code(&dec, LZMA_FINISH);
        if (ret == LZMA_STREAM_END)
            break;
        if (ret != LZMA_OK) {
            fprintf(stderr, "Decompression failed: %d\n", ret);
            free(compressed_buf);
            free(decompressed_buf);
            lzma_end(&dec);
            return 1;
        }
    }

    size_t decompressed_size = randbuf_len - dec.avail_out;
    lzma_end(&dec);

    /* Verify */
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
    printf("liblzma version: %s\n", LZMA_VERSION_STRING);
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
