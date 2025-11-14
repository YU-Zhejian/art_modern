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

#pragma once

#include "libam_support/Constants.hh"

#include <fmt/format.h>

#include <cstddef> // std::size_t
#include <string>

namespace labw::art_modern {
std::string format_with_commas(std::size_t number);

template <typename T> std::string to_si(T number, int precision = 2, T base = 1024)
{
    std::size_t unit_index = 0;
    std::string is_neg;
    auto size_in_unit = static_cast<double>(number);
    if (size_in_unit < 0) {
        is_neg = "-";
        size_in_unit = -size_in_unit;
    }

    while (size_in_unit >= base && unit_index < SI_UNITS_LENGTH) {
        size_in_unit /= base;
        ++unit_index;
    }
    // Named-arg variant (explicitly binds precision)
    return fmt::format("{0}{1:.{2}f}{3}", is_neg, size_in_unit, precision, SI_UNITS[unit_index]);
}
} // namespace labw::art_modern
