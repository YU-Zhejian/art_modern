#pragma once

#include "libam_support/Dtypes.h"

#include <cstdlib>
#include <tuple>

namespace labw::art_modern {

std::size_t num_base_positive_pe(
    am_read_len_t read_len_1, am_read_len_t read_len_2, am_readnum_t num_pos_reads, am_readnum_t num_neg_reads);
std::size_t num_base_negative_pe(
    am_read_len_t read_len_1, am_read_len_t read_len_2, am_readnum_t num_pos_reads, am_readnum_t num_neg_reads);
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads_template(
    double cov_pos, double cov_neg, am_readnum_t num_reads_to_reduce);
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads_se(
    std::size_t contig_size, am_read_len_t read_len_1, double cov_pos, double cov_neg);
/**
 * TODO: Slow method that used a gradient-descent-like approach. Optimize it.
 */
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads_pe(
    std::size_t contig_size, am_read_len_t read_len_1, am_read_len_t read_len_2, double cov_pos, double cov_neg);
} // namespace labw::art_modern
