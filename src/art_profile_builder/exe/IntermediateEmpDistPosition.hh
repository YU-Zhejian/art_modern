#pragma once

#include "libam_support/Constants.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace labw::art_modern {
class IntermediateEmpDistPosition {
public:
    constexpr static std::size_t ALL_IDX = 0;
    constexpr static std::size_t A_IDX = 1;
    constexpr static std::size_t C_IDX = 2;
    constexpr static std::size_t G_IDX = 3;
    constexpr static std::size_t T_IDX = 4;
    constexpr static std::size_t N_IDX = 5;
    constexpr static std::size_t BASE_IDX[] = { ALL_IDX, A_IDX, C_IDX, G_IDX, T_IDX, N_IDX };
    constexpr static char IDX_BASE[] = { '.', 'A', 'C', 'G', 'T', 'N' };
    constexpr static std::size_t NUM_BASES = 6; // ACGTN + all
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
    ~IntermediateEmpDistPosition() = default;
    DEFAULT_COPY(IntermediateEmpDistPosition)
    DELETE_MOVE(IntermediateEmpDistPosition)

    void add(char base, am_qual_t qual);

    /**
     * Call this before writing.
     */
    void accumulate();

    void add(IntermediateEmpDistPosition const& other);

    void write(std::ostream& oss, std::size_t pos_id, std::size_t base_idx) const;

private:
    std::vector<std::size_t> data_;
};

} // namespace labw::art_modern
