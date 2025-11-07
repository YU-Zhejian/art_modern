#include "art/lib/ArtParams.hh"

#include "art/lib/ArtConstants.hh"
#include "art/lib/Empdist.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/arithmetic_utils.hh"

#include <htslib/hts.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

namespace {

    static std::array<double, HIGHEST_QUAL> gen_err_prob_()
    {
        std::array<double, HIGHEST_QUAL> tmp_err_prob {};
        for (int i = 0; i < HIGHEST_QUAL; i++) {
            tmp_err_prob[i] = std::pow(10, -i / 10.0);
        }
        return tmp_err_prob;
    }
} // namespace
ArtParams::ArtParams(const labw::art_modern::SIMULATION_MODE art_simulation_mode,
    const labw::art_modern::ART_LIB_CONST_MODE art_lib_const_mode, const bool sep_flag, std::string&& id,
    const int max_n, const am_read_len_t read_len_1, const am_read_len_t read_len_2, const double pe_frag_dist_mean,
    const double pe_frag_dist_std_dev, std::vector<double>&& per_base_ins_rate_1,
    std::vector<double>&& per_base_del_rate_1, std::vector<double>&& per_base_ins_rate_2,
    std::vector<double>&& per_base_del_rate_2, labw::art_modern::Empdist&& qdist,
    const std::size_t job_pool_reporting_interval_seconds,
    const std::size_t art_job_executor_reporting_interval_seconds)
    : art_simulation_mode(art_simulation_mode)
    , art_lib_const_mode(art_lib_const_mode)
    , sep_flag(sep_flag)
    , id(std::move(id))
    , max_n(max_n)
    , read_len_1(read_len_1)
    , read_len_2(read_len_2)
    , pe_frag_dist_mean(pe_frag_dist_mean)
    , pe_frag_dist_std_dev(pe_frag_dist_std_dev)
    , per_base_ins_rate_1(std::move(per_base_ins_rate_1))
    , per_base_del_rate_1(std::move(per_base_del_rate_1))
    , per_base_ins_rate_2(std::move(per_base_ins_rate_2))
    , per_base_del_rate_2(std::move(per_base_del_rate_2))
    , qdist(std::move(qdist))
    , job_pool_reporting_interval_seconds(job_pool_reporting_interval_seconds)
    , art_job_executor_reporting_interval_seconds(art_job_executor_reporting_interval_seconds)
    , err_prob(gen_err_prob_())
    , pe_dist_mean_minus_2_std(static_cast<hts_pos_t>(pe_frag_dist_mean - 2 * pe_frag_dist_std_dev))
    , read_len_max(am_max(read_len_1, read_len_2))
{
}

} // namespace labw::art_modern
