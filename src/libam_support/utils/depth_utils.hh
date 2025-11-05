#pragma once

#include "libam_support/Dtypes.h"

#include <cstdlib>
#include <tuple>

namespace labw::art_modern {
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads_template(
    double cov_pos, double cov_neg, am_readnum_t num_reads_to_reduce);
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads(std::size_t contig_size, am_read_len_t read_len_1,
    am_read_len_t read_len_2, double cov_pos, double cov_neg, am_readnum_t num_reads_to_reduce);
} // namespace labw::art_modern
