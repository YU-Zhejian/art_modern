#pragma once

#define PARALLEL_DISABLE (-1)
#define PARALLEL_ALL 0
#define ALN_GAP '-'

#define ART_VERSION "2.5.8"
#define ART_MODERN_VERSION "1.0.0"
#define MAPQ_MAX ((1 << 8) - 1) // 255

enum CIGAR {
    ALN_MATCH = 'M',
    INSERTION = 'I',
    DELETION = 'D',
    SKIP = 'N',
    SCLIP = 'S',
    HCLIP = 'H',
    PADDING = 'P',
    SEQ_MATCH = '=',
    SEQ_MISMATCH = 'X'

    // M 0 alignment match (can be a sequence match or mismatch) yes yes
    // I 1 insertion to the reference yes no
    // D 2 deletion from the reference no yes
    // N 3 skipped region from the reference no yes
    // S 4 soft clipping (clipped sequences present in SEQ) yes no
    // H 5 hard clipping (clipped sequences NOT present in SEQ) no no
    // P 6 padding (silent deletion from padded reference) no no
    //= 7 sequence match yes yes
    // X 8 sequence mismatch yes yes
};
