#pragma once

#include "art_modern_config.h"

#define PARALLEL_DISABLE (-1)
#define PARALLEL_ALL 0
#define ALN_GAP '-'

#define ART_VERSION "2.5.8"
#define MAPQ_MAX ((1 << 8) - 1) // 255
#define PHRED_OFFSET 33

enum class SIMULATION_MODE { WGS, TRANS, TEMPLATE };
#define SIMULATION_MODE_WGS "wgs"
#define SIMULATION_MODE_TRANS "trans"
#define SIMULATION_MODE_TEMPLATE "template"

enum class INPUT_FILE_TYPE { FASTA, PBSIM3_TEMPLATE, AUTO };
#define INPUT_FILE_TYPE_FASTA "fasta"
#define INPUT_FILE_TYPE_PBSIM3_TEMPLATE "pbsim3_template"
#define INPUT_FILE_TYPE_AUTO "auto"
enum class INPUT_FILE_PARSER { MEMORY, HTSLIB, STREAM, AUTO };
#define INPUT_FILE_PARSER_MEMORY "memory"
#define INPUT_FILE_PARSER_HTSLIB "htslib"
#define INPUT_FILE_PARSER_STREAM "stream"
#define INPUT_FILE_PARSER_AUTO "auto"

#undef CEU_CM_IS_DEBUG // FIXME: CMake seems not working!

#define TID_FOR_UNMAPPED (-1)
