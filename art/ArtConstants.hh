#pragma once

const int HIGHEST_QUAL = 80;

const int MAX_DIST_NUMBER = 1000000;
const double DEFAULT_INS_RATE_1 = 0.00009;
const double DEFAULT_DEL_RATE_1 = 0.00011;
const double DEFAULT_INS_RATE_2 = 0.00015;
const double DEFAULT_DEL_RATE_2 = 0.00023;
const int DEFAULT_MAX_INDEL = (-1);
const int DEFAULT_MAX_NUM_N = (-1);
const int DEFAULT_BATCH_SIZE = (1 << 14);

enum class ART_LIB_CONST_MODE { SE, PE, MP };
const char ART_LIB_CONST_MODE_SE[] = "se";
const char ART_LIB_CONST_MODE_PE[] = "pe";
const char ART_LIB_CONST_MODE_MP[] = "mp";
const char ART_PROGRAM_NAME[] = "art_modern";
