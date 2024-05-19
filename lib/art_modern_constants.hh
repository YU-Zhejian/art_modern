#pragma once
#include <ceu_check/ceu_check_os_macro.h>

// Following constants copied from HTSLib <sam.h>.
/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED 1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR 2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP 4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP 8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE 16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE 32
/*! @abstract this is read1 */
#define BAM_FREAD1 64
/*! @abstract this is read2 */
#define BAM_FREAD2 128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY 256
/*! @abstract QC failure */
#define BAM_FQCFAIL 512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP 1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

#define PARALLEL_DISABLE (-1)
#define PARALLEL_ALL 0
#define ALN_GAP '-'

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

#ifndef NULL_DEVICE
#ifdef CEU_ON_WINDOWS_64
#define NULL_DEVICE "NUL"
#else
#define NULL_DEVICE "/dev/null"
#endif
#endif
