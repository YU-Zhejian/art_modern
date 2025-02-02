#pragma once

#include "libam/Dtypes.hh"

#include <htslib/hts.h>

#include <cstddef>

namespace labw::art_modern {
std::size_t hts_tell(htsFile* fp);
} // namespace labw::art_modern
