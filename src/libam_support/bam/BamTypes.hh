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

#include "htslib/sam.h"

#include <memory>

namespace labw::art_modern {

struct BamDestroyer {
    void operator()(bam1_t* b) const { bam_destroy1(b); }
};

/**
 * Unique pointer of bam1_t data types.
 */
using bam1_t_uptr = std::unique_ptr<bam1_t, BamDestroyer>;
} // namespace labw::art_modern
