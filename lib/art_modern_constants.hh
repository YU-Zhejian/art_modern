#pragma once

#include "art_modern_config.h"

#define PARALLEL_DISABLE (-1)
#define PARALLEL_ALL 0
#define ALN_GAP '-'

#define ART_VERSION "2.5.8"
#define MAPQ_MAX ((1 << 8) - 1) // 255
#define PHRED_OFFSET 33
enum class SIMULATION_MODE { WGS, TRANS, TEMPLATE };
enum class INPUT_FILE_TYPE { FASTA, PBSIM3_TEMPLATE };
;
enum class INPUT_FILE_PARSER { MEMORY, HTSLIB, STREAM };
