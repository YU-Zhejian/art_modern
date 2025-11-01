/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "libam_support/utils/si_utils.hh"

#include <fmt/format.h>

#include <cstddef>
#include <cstdint>
#include <string>

namespace labw::art_modern {
std::string format_with_commas(const std::size_t number)
{
    std::string num_str = std::to_string(number);
    int64_t insertPosition = static_cast<int64_t>(num_str.length()) - 3;

    while (insertPosition > 0) {
        num_str.insert(insertPosition, ",");
        insertPosition -= 3;
    }

    return num_str;
}

} // namespace labw::art_modern
