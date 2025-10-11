#pragma once

#include <cstdlib>
#include <string>

namespace labw::art_modern {
/**
 * Art Profile Builder Configuration
 */
struct APBConfig {
    const std::string input_file_path;
    const std::size_t read_length;
    const std::size_t num_threads;
    const std::size_t num_io_threads;
    const bool is_pe;
    const bool is_ob;
    const std::string output_1_file_path;
    const std::string output_2_file_path;
};

} // namespace labw::art_modern
