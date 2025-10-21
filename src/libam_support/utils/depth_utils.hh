#pragma once

#include "libam_support/Dtypes.hh"

#include <cstdlib>
#include <tuple>

namespace labw::art_modern {
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads(
    std::size_t contig_size, int read_len, double cov_pos, double cov_neg, am_readnum_t num_reads_to_reduce);
} // namespace labw::art_modern
