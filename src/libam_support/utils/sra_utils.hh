#pragma once

#include <string>

namespace labw::art_modern {
constexpr char SRA_SIGNATURE[8] = { 'N', 'C', 'B', 'I', '.', 's', 'r', 'a' };
bool detect_sra(const std::string& file_path);
void assert_is_sra(const std::string& file_path);
} // namespace labw::art_modern
