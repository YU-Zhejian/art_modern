/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "art_profile_builder/lib/IntermediateEmpDistPosition.hh"

#include "libam_support/Dtypes.h"
#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/sam.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace labw::art_modern {
/**
 * An intermediate empirical distribution for a set of reads.
 */
class IntermediateEmpDist {
public:
    explicit IntermediateEmpDist(std::size_t read_length);

    DEFAULT_COPY_CONSTRUCTOR(IntermediateEmpDist)
    DELETE_COPY_ASSIGNMENT(IntermediateEmpDist)
    DELETE_MOVE(IntermediateEmpDist)
    DEFAULT_DESTRUCTOR(IntermediateEmpDist)

    bool parse_read(const bam1_t* b);
    bool parse_read(const std::string& seq, const std::vector<am_qual_t>& qual);
    bool parse_read_ngs_fragment(const std::string& seq, const std::string& qual_str);

    void accumulate();

    void add(const IntermediateEmpDist& other);

    void write(std::ostream& oss, bool is_ob) const;
    [[nodiscard]] std::size_t get_total_reads() const;
    [[nodiscard]] std::size_t get_total_bases() const;

private:
    std::vector<IntermediateEmpDistPosition> positions_;
    const std::size_t read_length_;
    std::size_t total_reads_ { 0 };
    std::size_t total_bases_ { 0 };
};
} // namespace labw::art_modern
