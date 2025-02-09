#pragma once

#include "libam_support/Dtypes.hh"

#include <string>

namespace labw::art_modern {
void validate_input_filename(const std::string& input_file_path, const std::string& arg_name);

am_filelen_t get_file_size(const std::string& file_path) noexcept;

void ensure_directory_exists(const std::string& dir_path);
void prepare_writer(const std::string& output_file_path);

} // namespace labw::art_modern
