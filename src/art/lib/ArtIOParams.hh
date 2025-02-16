#pragma once
#include "libam_support/Constants.hh"
#include "libam_support/ds/CoverageInfo.hh"

#include <boost/program_options/variables_map.hpp>

#include <string>
#include <vector>

namespace labw::art_modern {

struct ArtIOParams {
    const std::string input_file_name;
    const INPUT_FILE_TYPE art_input_file_type;
    const INPUT_FILE_PARSER art_input_file_parser;
    const CoverageInfo coverage_info;
    const int parallel;
    const int batch_size;
    const boost::program_options::variables_map vm;
    const std::vector<std::string> args;
};
} // namespace labw::art_modern
