#pragma once
#include "ArtConstants.hh"
#include "Empdist.hh"
#include "art_modern_constants.hh"
#include "fasta/CoverageInfo.hh"
#include "out/BaseReadOutput.hh"
#include <array>
#include <boost/program_options.hpp>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace labw {
namespace art_modern {

    struct ArtParams {
        SIMULATION_MODE art_simulation_mode;
        ART_LIB_CONST_MODE art_lib_const_mode;
        std::string input_file_name;
        INPUT_FILE_TYPE art_input_file_type;
        INPUT_FILE_PARSER art_input_file_parser;
        int parallel;
        bool sep_flag;
        std::string id;
        CoverageInfo coverage_info;
        int read_len;
        double pe_frag_dist_mean;
        double pe_frag_dist_std_dev;
        std::vector<double> per_base_ins_rate_1;
        std::vector<double> per_base_del_rate_1;
        std::vector<double> per_base_ins_rate_2;
        std::vector<double> per_base_del_rate_2;
        std::array<double, HIGHEST_QUAL> err_prob;
        hts_pos_t pe_dist_mean_minus_2_std;
        Empdist qdist;
        std::shared_ptr<BaseFastaFetch> fasta_fetch;
        std::shared_ptr<BaseReadOutput> out_dispatcher;
    };

} // namespace art_modern
} // namespace labw
