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

#include "libam_support/ds/CoverageInfo.hh"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <istream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
CoverageInfo::CoverageInfo(const double static_coverage)
    : static_coverage_positive_(static_coverage / 2)
    , static_coverage_negative_(static_coverage / 2)
{
}

CoverageInfo::CoverageInfo(coverage_map&& coverage_positive, coverage_map&& coverage_negative)
    : static_coverage_positive_(0)
    , static_coverage_negative_(0)
    , coverage_positive_(std::move(coverage_positive))
    , coverage_negative_(std::move(coverage_negative))
{
}

double CoverageInfo::coverage_positive(const std::string& contig_name) const
{
    if (coverage_positive_.empty()) {
        return static_coverage_positive_;
    }

    if (const auto& find_result = coverage_positive_.find(contig_name); find_result != coverage_positive_.end()) {
        return find_result->second;
    }
    return 0.0;
}

double CoverageInfo::coverage_negative(const std::string& contig_name) const
{
    if (coverage_negative_.empty()) {
        return static_coverage_negative_;
    }

    if (const auto& find_result = coverage_negative_.find(contig_name); find_result != coverage_negative_.end()) {
        return find_result->second;
    }
    return 0.0;
}

CoverageInfo::CoverageInfo(std::istream& istream)
: static_coverage_positive_(0)
    , static_coverage_negative_(0)
{
    std::string buff;
    std::vector<std::string> tokens;
    while (!istream.eof()) {
        std::getline(istream, buff);
        if (buff.empty() || buff.front() == '#') {
            continue;
        }
        split(tokens, buff, boost::is_any_of("\t"));
        if (tokens.size() == 3) {
            coverage_positive_.try_emplace(tokens.at(0), std::stod(tokens.at(1)));
            coverage_negative_.try_emplace(tokens.at(0), std::stod(tokens.at(2)));
        } else if (tokens.size() == 2) {
            coverage_positive_.try_emplace(tokens.at(0), std::stod(tokens.at(1)) / 2);
            coverage_negative_.try_emplace(tokens.at(0), std::stod(tokens.at(1)) / 2);
        }
    }
}

CoverageInfo CoverageInfo::div(const int num_parts) const
{
    if (coverage_positive_.empty()) {
        return CoverageInfo(static_coverage_positive_ / num_parts, static_coverage_negative_ / num_parts);
    }

    coverage_map coverage_positive_new;
    for (auto const& [contig_name, contig_cov] : coverage_positive_) {
        coverage_positive_new.try_emplace(contig_name, contig_cov / num_parts);
    }
    coverage_map coverage_negative_new;
    for (auto const& [contig_name, contig_cov] : coverage_positive_) {
        coverage_negative_new.try_emplace(contig_name, contig_cov / num_parts);
    }
    return { std::move(coverage_positive_new), std::move(coverage_negative_new) };
}

CoverageInfo::CoverageInfo(const double static_coverage_positive, const double static_coverage_negative)
    : static_coverage_positive_(static_coverage_positive)
    , static_coverage_negative_(static_coverage_negative)
{
}


} // namespace labw::art_modern
