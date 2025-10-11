#include "art_profile_builder/lib/IntermediateEmpDistPosition.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace labw::art_modern {

IntermediateEmpDistPosition::IntermediateEmpDistPosition()
{
    data_.resize(NUM_BASES * WIDTH, 0);
    std::fill(data_.begin(), data_.end(), 0);
}

void IntermediateEmpDistPosition::add(const char base, const am_qual_t qual)
{
    std::size_t const base_idx = BASE_ASCII_TO_IDX[static_cast<unsigned char>(base)];
    if (qual > MAX_QUAL || qual < MIN_QUAL) {
        // Ignored
        return;
    }
    auto const qual_idx = static_cast<unsigned char>(qual - MIN_QUAL);
    data_[base_idx * WIDTH + qual_idx]++;
    data_[ALL_IDX * WIDTH + qual_idx]++;
}

void IntermediateEmpDistPosition::accumulate()
{
    for (std::size_t base_idx = 0; base_idx < NUM_BASES; ++base_idx) {
        for (std::size_t qual_idx = MIN_QUAL + 1; qual_idx <= MAX_QUAL; ++qual_idx) {
            data_[base_idx * WIDTH + qual_idx] += data_[base_idx * WIDTH + qual_idx - 1];
        }
    }
}

void IntermediateEmpDistPosition::add(IntermediateEmpDistPosition const& other)
{
    for (std::size_t i = 0; i < NUM_BASES * WIDTH; ++i) {
        data_[i] += other.data_[i];
    }
}

void IntermediateEmpDistPosition::write(
    std::ostream& oss, const std::size_t pos_id, const std::size_t base_idx, const bool is_ob) const
{
    const auto leading_char = IDX_BASE[base_idx];
    std::vector<std::size_t> output_qual_idx;
    std::size_t starting_pos = MAX_QUAL + 1;
    for (std::size_t qual_idx = MIN_QUAL; qual_idx <= MAX_QUAL; ++qual_idx) {
        if (data_[base_idx * WIDTH + qual_idx] != 0) {
            starting_pos = qual_idx;
            break;
        }
    }
    if (starting_pos != MAX_QUAL + 1) {
        output_qual_idx.emplace_back(starting_pos);
        for (std::size_t qual_idx = starting_pos + 1; qual_idx <= MAX_QUAL; ++qual_idx) {
            if (data_[base_idx * WIDTH + qual_idx] - data_[base_idx * WIDTH + qual_idx - 1] != 0) {
                output_qual_idx.emplace_back(qual_idx);
            }
        }
    }

    oss << leading_char << '\t' << pos_id << '\t';
    for (const auto qual_idx : output_qual_idx) {
        oss << qual_idx + (is_ob ? 1 : 0) << '\t'; // Here, quality '!' is denoted as 1
    }
    oss << '\n';

    oss << leading_char << '\t' << pos_id << '\t';
    for (const auto qual_idx : output_qual_idx) {
        oss << data_[base_idx * WIDTH + qual_idx] << '\t';
    }
    oss << '\n';
}
} // namespace labw::art_modern
