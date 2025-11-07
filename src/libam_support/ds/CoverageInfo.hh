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

#include "libam_support/utils/class_macros_utils.hh"

#include <istream>
#include <string>
#include <unordered_map>

namespace labw::art_modern {

class CoverageInfo {

public:
    using coverage_map = std::unordered_map<std::string, double>;

    explicit CoverageInfo(double static_coverage);
    explicit CoverageInfo(double static_coverage_positive, double static_coverage_negative);
    CoverageInfo(coverage_map&& coverage_positive, coverage_map&& coverage_negative);
    explicit CoverageInfo(std::istream& istream);

    DEFAULT_COPY(CoverageInfo)
    DEFAULT_MOVE(CoverageInfo)
    DEFAULT_DESTRUCTOR(CoverageInfo)

    [[nodiscard]] double coverage_positive(const std::string& contig_name) const;
    [[nodiscard]] double coverage_negative(const std::string& contig_name) const;
    [[nodiscard]] CoverageInfo div(std::size_t num_parts) const;

private:
    double static_coverage_positive_;
    double static_coverage_negative_;
    coverage_map coverage_positive_;
    coverage_map coverage_negative_;
};

} // namespace labw::art_modern
