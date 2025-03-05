#pragma once

#include <htslib/hts.h>

#include <cstddef>

namespace labw::art_modern {
/**
 * @brief Know where we are.
 * @param fp HTS file pointer
 * @return Where we are.
 */
std::size_t hts_tell(htsFile* fp);
} // namespace labw::art_modern
