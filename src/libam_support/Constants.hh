/**
 * @brief Constants used acceoss the library.
 *
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#pragma once
#include "libam_support/Dtypes.hh"

#include <cstddef>
#include <cstdint>
#include <string>

namespace labw::art_modern {

/** Parallelization diabled. */
constexpr int PARALLEL_DISABLE = -1;

/** Parallelization with all cores. */
constexpr int PARALLEL_ALL = 0;

/**
 * Alignment gap character.
 */
constexpr char ALN_GAP = '-';

/**
 * Maximum phred-based quality score.
 * 10^{-40} error rate, which is quite low.
 */
constexpr am_qual_t MAX_QUAL = 40;
/** Minimum phred-based quality score. */
constexpr am_qual_t MIN_QUAL = 0;

/** Size of 1K under SI units. */
constexpr int K_SIZE = 1U << 10U;
/** Size of 1M under SI units. */
constexpr int M_SIZE = 1U << 20U;
/** Size of 1G under SI units. */
constexpr int G_SIZE = 1U << 30U;

/**
 * Version of ART this program is based on.
 */
constexpr char ART_VERSION[] = "2.5.8";

/** Maximum mapping quality */
constexpr int MAPQ_MAX = (1 << 8) - 1; // 255

constexpr am_qual_t PHRED_OFFSET = 33;

enum class SIMULATION_MODE : std::int8_t { WGS, TRANS, TEMPLATE };

constexpr char SIMULATION_MODE_WGS[] = "wgs";

constexpr char SIMULATION_MODE_TRANS[] = "trans";

constexpr char SIMULATION_MODE_TEMPLATE[] = "template";

enum class INPUT_FILE_TYPE : std::int8_t { FASTA, PBSIM3_TRANSCRIPTS };

constexpr char INPUT_FILE_TYPE_FASTA[] = "fasta";

constexpr char INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS[] = "pbsim3_transcripts";

constexpr char INPUT_FILE_TYPE_AUTO[] = "auto";

enum class INPUT_FILE_PARSER : std::int8_t { MEMORY, HTSLIB, STREAM };

constexpr char INPUT_FILE_PARSER_MEMORY[] = "memory";

constexpr char INPUT_FILE_PARSER_HTSLIB[] = "htslib";

constexpr char INPUT_FILE_PARSER_STREAM[] = "stream";

constexpr char INPUT_FILE_PARSER_AUTO[] = "auto";

constexpr int MPI_MAIN_RANK = 0;
const std::string MPI_MAIN_RANK_STR = "0";

constexpr int MPI_UNAVAILABLE_RANK = -1;

constexpr int TID_FOR_UNMAPPED = -1;

/**
 * BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF; Consume both
 */
constexpr am_cigar_type_t CONSUME_QUERY_AND_REFERENCE = 0b11;

/**
 * BAM_CHARD_CLIP, BAM_CPAD, BAM_CBACK; Consume neither
 */
constexpr am_cigar_type_t CONSUME_NEITHER_QUERY_NOR_REFERENCE = 0b00;

/**
 *  BAM_CINS, BAM_CSOFT_CLIP; Consume query
 */
constexpr am_cigar_type_t CONSUME_QUERY = 0b01;

/**
 *  BAM_CDEL, BAM_CREF_SKIP; Consume reference
 */
constexpr am_cigar_type_t CONSUME_REFERENCE = 0b10;

/**
 * SI units for printing
 */
constexpr const char* SI_UNITS[] = { "", "K", "M", "G", "T", "P", "E", "Z", "Y" };
constexpr std::size_t SI_UNITS_LENGTH = 8;

/** ART nucleotides */
constexpr char ART_ACGT[] = "ACGT";
/** ART nucleotides as string */
const std::string ART_ACGT_STR = "ACGT";
} // namespace labw::art_modern
