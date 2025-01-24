#pragma once
#include "art/ArtConstants.hh"
#include "art/Empdist.hh"

#include "libam/Constants.hh"
#include "libam/ds/CoverageInfo.hh"

#include <absl/base/attributes.h>

#include <boost/program_options/variables_map.hpp>

#include <htslib/hts.h>

#include <array>
#include <string>
#include <vector>

namespace labw::art_modern {

struct ArtParams {
    const SIMULATION_MODE art_simulation_mode;
    const ART_LIB_CONST_MODE art_lib_const_mode;
    const std::string input_file_name;
    const INPUT_FILE_TYPE art_input_file_type;
    const INPUT_FILE_PARSER art_input_file_parser;
    const int parallel;
    const bool sep_flag;
    const std::string id;
    const CoverageInfo coverage_info;
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
    const int batch_size;
    const boost::program_options::variables_map vm;
    const std::vector<std::string> args;
} ABSL_ATTRIBUTE_PACKED;

} // namespace labw::art_modern
