#pragma once
#include "art_modern_dtypes.hh"
#include <string>
namespace labw::art_modern {

constexpr am_qual_t HIGHEST_QUAL = 80;

constexpr am_qual_dist_t MAX_DIST_NUMBER = 1000000;
constexpr double DEFAULT_INS_RATE_1 = 0.00009;
constexpr double DEFAULT_DEL_RATE_1 = 0.00011;
constexpr double DEFAULT_INS_RATE_2 = 0.00015;
constexpr double DEFAULT_DEL_RATE_2 = 0.00023;
constexpr int DEFAULT_MAX_INDEL = -1;
constexpr int DEFAULT_MAX_N = 0;
constexpr int DEFAULT_BATCH_SIZE = 1 << 14;

enum class ART_LIB_CONST_MODE { SE, PE, MP };
constexpr char ART_LIB_CONST_MODE_SE[] = "se";
constexpr char ART_LIB_CONST_MODE_PE[] = "pe";
constexpr char ART_LIB_CONST_MODE_MP[] = "mp";
constexpr char ART_PROGRAM_NAME[] = "art_modern";
constexpr char ART_ACGT[] = "ACGT";
const std::string ART_ACGT_STR = "ACGT";

// If 20% attempts failed
constexpr double MAX_TRIAL_RATIO_BEFORE_FAIL = 0.2;

}