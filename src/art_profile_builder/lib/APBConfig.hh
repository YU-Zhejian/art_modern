#pragma once

#include "libam_support/Dtypes.h"

#include <htslib/hts.h>

#include <cstdlib>
#include <string>

namespace labw::art_modern {
/**
 * Art Profile Builder Configuration
 */
class APBConfig {
public:
    const std::string input_file_path;
    const am_read_len_t read_length_1;
    const am_read_len_t read_length_2;
    const std::size_t num_threads;
    const std::size_t num_io_threads;
    const bool is_pe;
    const bool is_ob;
    const std::string output_1_file_path;
    const std::string output_2_file_path;
    const htsExactFormat input_file_format;
};

} // namespace labw::art_modern
