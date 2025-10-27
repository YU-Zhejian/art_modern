
#include <libdeflate.h>
#include <zlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getdelim.h"

const size_t NUM_CLASSES = 16 * 16; // 256 classes for 2 hex digits

int main(int argc, char** argv)
{
    // The main process will read line from /dev/stdin.
    // It will then classify the line using the first 2 hex digit of its CRC32 checksum.
    // Finally, it will print the classification result to file named [basename].[classification].
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <basename> < /dev/stdin > /dev/stdout\n", argv[0]);
        return EXIT_FAILURE;
    }
    const char* basename = argv[1];

    gzFile* output_files = (gzFile*)calloc(NUM_CLASSES, sizeof(gzFile));
    if (output_files == NULL) {
        perror("calloc");
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < NUM_CLASSES; i++) {
        char filename[256];
        snprintf(filename, sizeof(filename), "%s.%02zx.gz", basename, i);
        output_files[i] = gzopen(filename, "w");
        if (output_files[i] == NULL) {
            perror("gzopen");
            return EXIT_FAILURE;
        }
    }

    char* line = NULL;
    size_t len = 0;
    size_t n_lines = 0;
    while ((getline(&line, &len, stdin)) != -1) {
        uint32_t crc32 = libdeflate_crc32(0, line, len);
        size_t class = (crc32 >> 4) & 0xFF;
        gzprintf(output_files[class], "%s", line);
        n_lines++;
        if (n_lines % 1000000 == 0) {
            fprintf(stderr, "Processed %zu lines...\n", n_lines);
            fflush(stderr);
        }
    }
    free(line);
    for (size_t i = 0; i < NUM_CLASSES; i++) {
        gzflush(output_files[i], Z_SYNC_FLUSH);
        gzclose(output_files[i]);
    }
    free(output_files);
    return EXIT_SUCCESS;
}
