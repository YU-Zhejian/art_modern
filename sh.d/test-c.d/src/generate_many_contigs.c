#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
    const size_t CONTIG_SIZE = 1024ULL;
    const size_t NUM_CONTIGS = 5ULL * 1024 * 1024 * 1024;
    const size_t LINE_LEN = CONTIG_SIZE;
    // Generate 1M random numbers for later use
    char line[LINE_LEN + 1];
    memset(line, 'A', LINE_LEN);
    line[LINE_LEN] = '\0';
    for (size_t i = 0; i < NUM_CONTIGS; i++) {
        printf(">contig%09lx\n%s\n", i, line);
    }
    return EXIT_SUCCESS;
}
