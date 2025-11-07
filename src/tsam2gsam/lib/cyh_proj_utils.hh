#pragma once

#include "libam_support/Dtypes.h"

#include <htslib/sam.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace labw::art_modern {

/**
 * Version of the current programs.
 */
constexpr char TSAM2GSAM_VERSION[] = "2025-11-06";

/**
 * Number of exons expected when initializing the list.
 */
constexpr int NUM_EXONS_EXPECTED = 10;

// GFFUtils BED are bed12 format with 12 fields.
constexpr int NUM_GFFUTILS_BED_FIELDS = 12;

/**
 * Indicating sam_read1 reaching end of read.
 */
constexpr int SAM_READ_EOF = -1;

/** Allowed number of bases in splice sites **/
constexpr int ALLOWED_ERROR_BASES = 5;

constexpr int N_THREADS_FOR_HTSLIB = 10;

/*!
 * Merge CIGARs using the following rules:
 *
 * - If the first or last CIGAR operation is an insertion, convert it to soft clipping.
 * - Merge adjacent CIGAR operations of the same type.
 * - Ignore CIGAR operations with zero length.
 *
 * @param g_cigar Established cigars.
 * @return
 */
std::vector<am_cigar_t> merge_cigars(const std::vector<am_cigar_t>& g_cigar);

/**
 * Convert BAM quality to C++ string.
 *
 * Some sucker may not provide any quality.
 * If such, generate a dumb quality string filled with !.
 *
 * @param aln The BAM alignment.
 * @return Converted quality string.
 */
std::string bam_qual_to_str(const bam1_t* aln);

/**
 * Clear pair-end mapping flag from BAM alignment.
 *
 * @param flag Mutable reference to BAM alignment flag.
 */
void clear_pe_flag(uint16_t& flag);

/**
 * Convert query sequence inside BAM record to string.
 *
 * @param aln BAM alignment.
 * @return The query sequence.
 */
std::string bam_seq_to_str(const bam1_t* aln);

/**
 * Set sequence inside BAM record. Note that this method can only be used if the new sequence is of equal length of the
 * old one.
 *
 * @param aln BAM alignment.
 * @param seq The query sequence.
 */
void bam_set_seq_eqlen(bam1_t* aln, const char* seq);

const char* get_sam_mode(const char* file_name, bool write);

void print_version(const std::string& prog_name);
} // namespace labw::art_modern
