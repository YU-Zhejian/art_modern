#pragma once
#include "art_modern_dtypes.hh"
#include <string>

namespace labw::art_modern {

constexpr int PARALLEL_DISABLE = -1;

constexpr int PARALLEL_ALL = 0;

constexpr char ALN_GAP = '-';

const static std::string ALN_GAP_STR = "-";

/**
 * 10^{-40} error rate, which is quite low.
 */
constexpr am_qual_t MAX_QUAL = 40;
constexpr am_qual_t MIN_QUAL = 0;
constexpr int K_SIZE = 1U << 10U;
constexpr int M_SIZE = 1U << 20U;
constexpr int G_SIZE = 1U << 30U;

constexpr char ART_VERSION[] = "2.5.8";

constexpr int MAPQ_MAX = ((1 << 8) - 1); // 255
constexpr char MAPQ_MAX_STR[] = "255"; // 255

constexpr int PHRED_OFFSET = 33;

enum class SIMULATION_MODE { WGS, TRANS, TEMPLATE };

constexpr char SIMULATION_MODE_WGS[] = "wgs";

constexpr char SIMULATION_MODE_TRANS[] = "trans";

constexpr char SIMULATION_MODE_TEMPLATE[] = "template";

enum class INPUT_FILE_TYPE { FASTA, PBSIM3_TRANSCRIPTS };

constexpr char INPUT_FILE_TYPE_FASTA[] = "fasta";

constexpr char INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS[] = "pbsim3_transcripts";

constexpr char INPUT_FILE_TYPE_AUTO[] = "auto";

enum class INPUT_FILE_PARSER { MEMORY, HTSLIB, STREAM };

constexpr char INPUT_FILE_PARSER_MEMORY[] = "memory";

constexpr char INPUT_FILE_PARSER_HTSLIB[] = "htslib";

constexpr char INPUT_FILE_PARSER_STREAM[] = "stream";

constexpr char INPUT_FILE_PARSER_AUTO[] = "auto";

constexpr int MPI_MAIN_RANK = 0;

constexpr int MPI_UNAVAILABLE_RANK = -1;

constexpr int TID_FOR_UNMAPPED = -1;

/**
 * BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF; Consume both
 */
constexpr int CONSUME_QUERY_AND_REFERENCE = 0b11;

/**
 * BAM_CHARD_CLIP, BAM_CPAD, BAM_CBACK; Consume neither
 */
constexpr int CONSUME_NEITHER_QUERY_NOR_REFERENCE = 0b00;

/**
 *  BAM_CINS, BAM_CSOFT_CLIP; Consume query
 */
constexpr int CONSUME_QUERY = 0b01;

/**
 *  BAM_CDEL, BAM_CREF_SKIP; Consume reference
 */
constexpr int CONSUME_REFERENCE = 0b10;
} // namespace labw::art_modern
