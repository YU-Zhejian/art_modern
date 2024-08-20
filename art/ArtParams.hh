#pragma once
#include "ArtConstants.hh"
#include "Empdist.hh"
#include "art_modern_constants.hh"
#include "out/OutputDispatcher.hh"
#include "fasta/CoverageInfo.hh"
#include <boost/program_options.hpp>
#include <map>
#include <string>
#include <vector>
#include <array>

namespace labw {
namespace art_modern {

    struct ArtParams {
        SIMULATION_MODE art_simulation_mode;
        ART_LIB_CONST_MODE art_lib_const_mode;
        INPUT_FILE_TYPE art_input_file_type;
        INPUT_FILE_PARSER art_input_file_parser;
        int parallel;
        bool sep_flag;
        std::string id;
        std::string seq_file;
        std::string fcov;
        int read_len;
        int pe_frag_dist_mean;
        double pe_frag_dist_std_dev;
        std::vector<double> per_base_ins_rate_1;
        std::vector<double> per_base_del_rate_1;
        std::vector<double> per_base_ins_rate_2;
        std::vector<double> per_base_del_rate_2;
        std::array<double, HIGHEST_QUAL> err_prob;
        hts_pos_t pe_dist_mean_minus_2_std;
        Empdist qdist_;
        CoverageInfo coverage_info;
        OutputDispatcherFactory out_dispatcher_factory_;
    };

} // namespace art_modern
} // namespace labw
