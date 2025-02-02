#pragma once

#include "htslib/sam.h"

#include <memory>

namespace labw::art_modern {

struct BamDestroyer {
    void operator()(bam1_t* b) const { bam_destroy1(b); }
};

/**
 * Unique pointer of bam1_t data types.
 */
using bam1_t_uptr = std::unique_ptr<bam1_t, BamDestroyer>;
} // namespace labw::art_modern