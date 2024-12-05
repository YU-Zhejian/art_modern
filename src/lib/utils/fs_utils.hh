#pragma once
#include <string>
namespace labw {
namespace art_modern {
    void validate_input_filename(const std::string& input_file_path, const std::string& arg_name);

    long get_file_size(const std::string& file_path) noexcept;

} // art_modern
} // labw
