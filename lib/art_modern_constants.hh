#pragma once

#include "art_modern_config.h"

const int PARALLEL_DISABLE = -1;
const int PARALLEL_ALL = 0;
const char ALN_GAP = '-';

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

enum class INPUT_FILE_TYPE { FASTA, PBSIM3_TEMPLATE };
const char INPUT_FILE_TYPE_FASTA[] = "fasta";
const char INPUT_FILE_TYPE_PBSIM3_TEMPLATE[] = "pbsim3_template";
const char INPUT_FILE_TYPE_AUTO[] = "auto";

enum class INPUT_FILE_PARSER { MEMORY, HTSLIB, STREAM };
const char INPUT_FILE_PARSER_MEMORY[] = "memory";
const char INPUT_FILE_PARSER_HTSLIB[] = "htslib";
const char INPUT_FILE_PARSER_STREAM[] = "stream";
const char INPUT_FILE_PARSER_AUTO[] = "auto";

#undef CEU_CM_IS_DEBUG // FIXME: CMake seems not working!

const int TID_FOR_UNMAPPED = (-1);
