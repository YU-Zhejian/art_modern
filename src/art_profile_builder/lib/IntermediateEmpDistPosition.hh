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

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstdlib>
#include <iostream>
#include <vector>

namespace labw::art_modern {
class IntermediateEmpDistPosition {
public:
    constexpr static std::size_t ALL_IDX = 0;
    constexpr static std::size_t A_IDX = 1;
    constexpr static std::size_t T_IDX = 2;
    constexpr static std::size_t G_IDX = 3;
    constexpr static std::size_t C_IDX = 4;
    constexpr static std::size_t N_IDX = 5;
    constexpr static std::size_t BASE_IDX[] = { ALL_IDX, A_IDX, T_IDX, G_IDX, C_IDX, N_IDX };
    constexpr static char IDX_BASE[] = { '.', 'A', 'T', 'G', 'C', 'N' };
    constexpr static std::size_t NUM_BASES = 6; // ACGTN + all
    constexpr static std::size_t WIDTH = (MAX_QUAL - MIN_QUAL + 1);
    /** ASCII to index
     *
     *  Generated using Python:
     *  @code
     *  for i in range(0, 256): print((chr(i).upper() if chr(i) in "ACGTacgt" else "N")+"_IDX", end=", ")
     *  @endcode
     */
    constexpr static std::size_t BASE_ASCII_TO_IDX[] = { N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, A_IDX, N_IDX, C_IDX, N_IDX, N_IDX, N_IDX, G_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, T_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, A_IDX, N_IDX, C_IDX, N_IDX, N_IDX, N_IDX, G_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, T_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX };

    IntermediateEmpDistPosition();
    DEFAULT_DESTRUCTOR(IntermediateEmpDistPosition)
    DEFAULT_COPY(IntermediateEmpDistPosition)
    DEFAULT_MOVE(IntermediateEmpDistPosition)

    void add(char base, am_qual_t qual);

    /**
     * Call this before writing.
     */
    void accumulate();

    void add(IntermediateEmpDistPosition const& other);

    void write(std::ostream& oss, std::size_t pos_id, std::size_t base_idx, bool is_ob) const;

private:
    std::vector<std::size_t> data_;
};

} // namespace labw::art_modern
