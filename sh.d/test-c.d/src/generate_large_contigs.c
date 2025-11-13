/*!
 * \file generate_large_contigs.c
 * \brief Generate large contigs for testing.
 * 
 * Note that the generated contigs are not biologically meaningful.
 * They are just sequences of 'A's to reach the desired size.
 * 
 * To generate random contigs, consider using:
 * - Use Intel MKL for random number generation. Generate uint32s.
 * - Map the uint32s to A/C/G/T using 2-bit decoding.
 * - Add telomeric Ns and centromeric Ns as needed.
 * - Use POSIX AIO for efficient writing.
 * 
 * TODO: Merge all C files into 1 directory.
 */


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
    return EXIT_SUCCESS;
}
