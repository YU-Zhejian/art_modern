/*!
 * \file generate_rand_large_contigs.c
 * \brief Generate random large contigs for testing.
 *
 * To generate random contigs, consider using:
 * - Use Intel MKL for random number generation. Generate uint32s.
 * - Map the uint32s to A/C/G/T using 2-bit decoding.
 * - Add telomeric Ns and centromeric Ns as needed.
 * - Use POSIX AIO for efficient writing.
 *
 * Speed: ~1.17 GiB/s.
 */

#include "twobit.h"

#include <mkl.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const size_t CONTIG_SIZE = 5ULL * 1024 * 1024 * 1024;
static const size_t CHUNK_SIZE = 4096ULL * 1024; // 4 MiB
static const size_t NUM_CONTIGS = 4;
static const size_t LINE_LEN = 64;

typedef struct cache {
    size_t n_rand_uint32s;
    uint32_t* rand_uint32s;
    char* nts;
} cache_t;

void generate_large_seq(twobit_decoder_t* twobit_decoder, VSLStreamStatePtr stream, cache_t* cache)
{
    const size_t n_rand_uint8s = cache->n_rand_uint32s << 2;
    const size_t n_bases = n_rand_uint8s << 2;

    viRngUniformBits32(VSL_RNG_METHOD_UNIFORM_STD, stream, (int)cache->n_rand_uint32s, cache->rand_uint32s);
    // Now cast uint32s to uint8s
    uint8_t* rand_uint8s = (uint8_t*)cache->rand_uint32s;
    // Now map the uint8s to nts using 2-bit decoding
    twobit_to_nt(twobit_decoder, rand_uint8s, cache->nts, 0, 0, n_rand_uint8s);
    // Now print the nts to stdout in lines of 64 chars
    for (size_t i = 0; i < n_bases / LINE_LEN; i++) {
        fwrite(cache->nts + i * LINE_LEN, sizeof(char), LINE_LEN, stdout);
        fputc('\n', stdout);
    }
}

int main(void)
{
    twobit_decoder_t* twobit_decoder = twobit_decoder_init();
    VSLStreamStatePtr stream = NULL;
    vslNewStream(&stream, VSL_BRNG_SFMT19937, 1);
    cache_t* cache = (cache_t*)malloc(sizeof(cache_t));
    cache->n_rand_uint32s = CHUNK_SIZE >> 4; // Since each uint32 gives 4 bytes, and each byte gives 4 nts.
    cache->rand_uint32s = (uint32_t*)mkl_malloc(cache->n_rand_uint32s * sizeof(uint32_t), 64);
    cache->nts = (char*)malloc((CHUNK_SIZE) * sizeof(char));

    for (size_t i = 0; i < NUM_CONTIGS; i++) {
        printf(">contig%ld\n", i);
        // Generate 32 GiB ATCGN sequence
        // That is, 64 chars per line for 512 M lines.
        for (size_t j = 0; j < CONTIG_SIZE / CHUNK_SIZE; j++) {
            generate_large_seq(twobit_decoder, stream, cache);
        }
    }

    mkl_free(cache->rand_uint32s);
    free(cache->nts);
    free(cache);
    twobit_decoder_free(twobit_decoder);
    vslDeleteStream(&stream);
    return EXIT_SUCCESS;
}
