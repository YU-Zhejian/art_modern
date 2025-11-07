#pragma once
#include "art_modern_config.h" // NOLINT: For CEU_CM_IS_DEBUG

#include "tsam2gsam/lib/Transcript.hh"

#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <ostream>

namespace labw::art_modern {
/**
 * @brief Convert transcript-aligned read to genome-aligned read.
 *
 * @code
 *
 *          Exon 0         Exon 1                    Exon 2
 *
 *  TALN   |<====|--------|====|--------------------|=============|
 *  GALN       |=|--------|====|--------------------|==>|
 *
 * @endcode
 *
 * @param t_aln Input transcript-aligned read.
 * @param g_aln Output genome-aligned read.
 * @param transcript The transcript.
 * @return true if there's no error; false otherwise
 */
#ifdef CEU_CM_IS_DEBUG
void convert_transcript_to_genome_alignment(
    bam1_t* t_aln, bam1_t* g_aln, const Transcript& transcript, std::ostream& cigar_trace);
#else
void convert_transcript_to_genome_alignment(bam1_t* t_aln, bam1_t* g_aln, const Transcript& transcript);
#endif
void populate_ghdr(sam_hdr_t* ghdr, faidx_t* faidx, int argc, char ** argv);
} // namespace labw::art_modern
