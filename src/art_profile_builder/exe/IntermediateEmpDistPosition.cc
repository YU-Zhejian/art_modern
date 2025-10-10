#include "art_profile_builder/exe/IntermediateEmpDistPosition.hh"

#include "libam_support/Constants.hh"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

namespace labw::art_modern {

IntermediateEmpDistPosition::IntermediateEmpDistPosition()
{
    data_.resize(NUM_BASES * (MAX_QUAL - MIN_QUAL + 1), 0);
    std::fill(data_.begin(), data_.end(), 0);
}

void IntermediateEmpDistPosition::add(const char base, const am_qual_t qual)
{
    std::size_t const base_idx = BASE_ASCII_TO_IDX[static_cast<unsigned char>(base)];
    if (qual > MAX_QUAL || qual < MIN_QUAL) {
        // Ignored
        return;
    }
    auto const qual_idx = static_cast<unsigned char>(qual);
    data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx]++;
    data_[ALL_IDX * (MAX_QUAL - MIN_QUAL + 1) + qual_idx]++;
}

void IntermediateEmpDistPosition::accumulate()
{
    for (std::size_t base_idx = 0; base_idx < NUM_BASES; ++base_idx) {
        for (std::size_t qual_idx = MIN_QUAL + 1; qual_idx <= MAX_QUAL; ++qual_idx) {
            data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx]
                += data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx - 1];
        }
    }
}

void IntermediateEmpDistPosition::add(IntermediateEmpDistPosition const& other)
{
    for (std::size_t i = 0; i < NUM_BASES * (MAX_QUAL - MIN_QUAL + 1); ++i) {
        data_[i] += other.data_[i];
    }
}

void IntermediateEmpDistPosition::write(std::ostream& oss, const std::size_t pos_id, const std::size_t base_idx) const
{
    char leading_char = IDX_BASE[base_idx];

    oss << leading_char << '\t' << pos_id << '\t';
    for (std::size_t qual_idx = MIN_QUAL; qual_idx <= MAX_QUAL; ++qual_idx) {
        if (data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx] != 0) {
            oss << qual_idx << '\t';
        }
    }
    oss << '\n';

    oss << leading_char << '\t' << pos_id << '\t';
    for (std::size_t qual_idx = MIN_QUAL; qual_idx <= MAX_QUAL; ++qual_idx) {
        if (data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx] != 0) {
            oss << data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx] << '\t';
        }
    }
    oss << '\n';
}
} // namespace labw::art_modern
