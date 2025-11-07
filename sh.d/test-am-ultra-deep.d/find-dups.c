
#include "dtypes.h"

#include <bzlib.h>
#include <xxhash.h>

#include <bloom/bloom.h>
#include <klib/khash.h>
// #include <netbsd/getdelim.h>
#include <hedley/hedley.h>
#include <rapidhash/rapidhash.h>

#define _POSIX_C_SOURCE 200809L // For getline, CLOCK_MONOTONIC

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <unistd.h>

KHASH_INIT(dup_set, kh_cstr_t, char, 0, kh_str_hash_func, kh_str_hash_equal)

const size_t BLOOM_FILTER_SIZE = 10ULL * 1024 * 1024 * 1024; // 10GiB
const size_t REPORT_INTERVAL = 1000000ULL;

hash_type xxh64_wrapper(const void* data, size_t len) { return XXH64(data, len, 0); }

hash_type rapidhash_wrapper(const void* data, size_t len) { return rapidhashMicro(data, len); }

int main(void)
{
    // Use a bloom filter to find duplicate lines from /dev/stdin.
    // Exit with code 0 if no duplicates are found.
    // Exit with code 1 if duplicates are found.

    // Initialize the bloom filter
    bloom_t bloom_filter = bloom_create(BLOOM_FILTER_SIZE);
    if (bloom_filter == NULL) {
        fprintf(stderr, "Failed to create bloom filter\n");
        return EXIT_FAILURE;
    }
    bloom_add_hash(bloom_filter, xxh64_wrapper);
    bloom_add_hash(bloom_filter, rapidhash_wrapper);

    // Initialize the khash set to store duplicates
    khash_t(dup_set)* dups = kh_init(dup_set);
    khiter_t k;
    if (dups == NULL) {
        fprintf(stderr, "Failed to create duplicate set\n");
        bloom_free(bloom_filter);
        return EXIT_FAILURE;
    }
    int32_t ret;

    // OK now read lines from stdin
    char* line = NULL;
    size_t bufflen = 0;
    ssize_t linelen;
    size_t linecount = 0;
    size_t bloom_errors = 0;
    bool fucked = false;
    struct timespec ts_start;
    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    while ((linelen = getline(&line, &bufflen, stdin)) != -1) {
        linecount++;
        // Check if the line is in the bloom filter
        if (true) { // bloom_test(bloom_filter, line, linelen)
            // Potential duplicate
            k = kh_get(dup_set, dups, line);
            if (k != kh_end(dups)) {
                // Confirmed duplicate
                fprintf(stderr, "Duplicate found (%zu): %s", linecount, line);
                fucked = true;
                goto cleanup;
            } else {
                bloom_errors++;
            }
        }
        // copy the line into the khash
        char* line_copy = strdup(line);
        kh_put(dup_set, dups, line_copy, &ret);
        bloom_add(bloom_filter, line, linelen);
        if (ret == -1) {
            // Failed to insert into khash
            fprintf(stderr, "Failed to insert into duplicate set\n");
            fucked = true;
            goto cleanup;
        }
        if (linecount % REPORT_INTERVAL == 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts_end);
            double elapsed = (ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec) / 1e9;
            double lps = REPORT_INTERVAL / elapsed;
            fprintf(stderr, "Read %zu lines; Bloom error rate: %f; Lines per second: %f\n", linecount,
                (double)bloom_errors / REPORT_INTERVAL, lps);
            bloom_errors = 0;
        }
    }
    fprintf(stderr, "Done\n");

cleanup:
    free(line);
    kh_destroy(dup_set, dups);
    bloom_free(bloom_filter);
    return fucked ? EXIT_FAILURE : EXIT_SUCCESS;
}
