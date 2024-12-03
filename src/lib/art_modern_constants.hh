#pragma once
#include <string>

const int PARALLEL_DISABLE = -1;
const int PARALLEL_ALL = 0;
const char ALN_GAP = '-';
const std::string ALN_GAP_STR = "-";

const int MAX_QUAL = 93;
const int MIN_QUAL = 0;

// !"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~

const char ART_VERSION[] = "2.5.8";
const int MAPQ_MAX = ((1 << 8) - 1); // 255
const int PHRED_OFFSET = 33;

enum class SIMULATION_MODE { WGS, TRANS, TEMPLATE };
const char SIMULATION_MODE_WGS[] = "wgs";
const char SIMULATION_MODE_TRANS[] = "trans";
const char SIMULATION_MODE_TEMPLATE[] = "template";

enum class INPUT_FILE_TYPE { FASTA, PBSIM3_TRANSCRIPTS };
const char INPUT_FILE_TYPE_FASTA[] = "fasta";
const char INPUT_FILE_TYPE_PBSIM3_TRANSCRIPTS[] = "pbsim3_transcripts";
const char INPUT_FILE_TYPE_AUTO[] = "auto";

enum class INPUT_FILE_PARSER { MEMORY, HTSLIB, STREAM };
const char INPUT_FILE_PARSER_MEMORY[] = "memory";
const char INPUT_FILE_PARSER_HTSLIB[] = "htslib";
const char INPUT_FILE_PARSER_STREAM[] = "stream";
const char INPUT_FILE_PARSER_AUTO[] = "auto";
const int MPI_MAIN_RANK = 0;
const int MPI_UNAVAILABLE_RANK = -1;

const int TID_FOR_UNMAPPED = -1;

/**
 * BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF; Consume both
 */
const int CONSUME_QUERY_AND_REFERENCE = 0b11;
/**
 * BAM_CHARD_CLIP, BAM_CPAD, BAM_CBACK; Consume neither
 */
const int CONSUME_NEITHER_QUERY_NOR_REFERENCE = 0b00;
/**
 *  BAM_CINS, BAM_CSOFT_CLIP; Consume query
 */
const int CONSUME_QUERY = 0b01;
/**
 *  BAM_CDEL, BAM_CREF_SKIP; Consume reference
 */
const int CONSUME_REFERENCE = 0b10;
