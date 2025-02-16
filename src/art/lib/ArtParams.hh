#pragma once
#include "art/lib/ArtConstants.hh"
#include "art/lib/Empdist.hh"

#include "libam_support/Constants.hh"

#include <htslib/hts.h>

#include <array>
#include <string>
#include <vector>

namespace labw::art_modern {

struct ArtParams {
    const SIMULATION_MODE art_simulation_mode;
    const ART_LIB_CONST_MODE art_lib_const_mode;
    const bool sep_flag;
    const std::string id;
    const int max_n;
    const int read_len;
    const double pe_frag_dist_mean;
    const double pe_frag_dist_std_dev;
    const std::vector<double> per_base_ins_rate_1;
    const std::vector<double> per_base_del_rate_1;
    const std::vector<double> per_base_ins_rate_2;
    const std::vector<double> per_base_del_rate_2;
    const std::array<double, HIGHEST_QUAL> err_prob;
    const hts_pos_t pe_dist_mean_minus_2_std;
    const Empdist qdist;
};

} // namespace labw::art_modern
