#include "depth_utils.hh"

#include "libam_support/Dtypes.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <tuple>

namespace labw::art_modern {
std::tuple<am_readnum_t, am_readnum_t> calculate_num_reads(const std::size_t contig_size, const int read_len,
    const double cov_pos, const double cov_neg, const am_readnum_t num_reads_to_reduce)
{
    const double cov_ratio = static_cast<double>(contig_size) / read_len;
    const auto joint_cov = cov_pos + cov_neg;

    const am_readnum_t num_c_pos_reads = std::lround(cov_pos * cov_ratio);
    const am_readnum_t num_c_neg_reads = std::lround(cov_neg * cov_ratio);
    const std::array<am_readnum_t, 3> num_pos_reads_candidates { num_c_pos_reads, num_c_pos_reads + 1,
        num_c_pos_reads - 1 };
    const std::array<am_readnum_t, 3> num_neg_reads_candidates { num_c_neg_reads, num_c_neg_reads - 1,
        num_c_neg_reads + 1 };
    std::size_t final_selected_pos_idx = 0;
    std::size_t final_selected_neg_idx = 0;
    am_readnum_t min_diff = std::numeric_limits<am_readnum_t>::max();
    const am_readnum_t num_total_bases = std::lround(joint_cov * static_cast<double>(contig_size));
    for (std::size_t i = 0; i < num_pos_reads_candidates.size(); ++i) {
        for (std::size_t j = 0; j < num_neg_reads_candidates.size(); ++j) {
            const am_readnum_t candidate_pos_reads = num_pos_reads_candidates[i];
            const am_readnum_t candidate_neg_reads = num_neg_reads_candidates[j];
            if (candidate_pos_reads < 0 || candidate_neg_reads < 0) {
                continue;
            }
            if (candidate_pos_reads % num_reads_to_reduce != 0 || candidate_neg_reads % num_reads_to_reduce != 0) {
                continue;
            }
            const am_readnum_t candidate_pos_bases = candidate_pos_reads * read_len;
            const am_readnum_t candidate_neg_bases = candidate_neg_reads * read_len;

            const am_readnum_t diff = (candidate_pos_bases + candidate_neg_bases) - num_total_bases;
            if (std::abs(diff) < std::abs(min_diff)) {
                min_diff = diff;
                final_selected_pos_idx = i;
                final_selected_neg_idx = j;
                if (diff == 0) {
                    break;
                }
            }
        }
    }
    const auto num_pos_reads = num_pos_reads_candidates[final_selected_pos_idx];
    const auto num_neg_reads = num_neg_reads_candidates[final_selected_neg_idx];
    return { num_pos_reads, num_neg_reads };
}
} // namespace labw::art_modern
