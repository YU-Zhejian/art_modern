#pragma once

#include <cstddef>
#include <string>

namespace labw::art_modern {
std::string format_with_commas(std::size_t number);
std::string to_si(double number, int precision = 2);
} // namespace labw::art_modern
