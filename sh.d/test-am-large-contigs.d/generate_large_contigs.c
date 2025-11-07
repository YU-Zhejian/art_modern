#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
    const size_t CONTIG_SIZE = 5ULL * 1024 * 1024 * 1024;
    const size_t NUM_CONTIGS = 4;
    const size_t LINE_LEN = 64;
    char line[LINE_LEN + 1];
    memset(line, 'A', LINE_LEN);
    line[LINE_LEN] = '\0';
    for (size_t i = 0; i < NUM_CONTIGS; i++) {
        printf(">contig%ld\n", i);
        // Generate 32 GiB ATCGN sequence
        // That is, 64 chars per line for 512 M lines.
        for (size_t j = 0; j < CONTIG_SIZE / LINE_LEN; j++) {
            printf("%s\n", line);
        }
    }
    return 0;
}
