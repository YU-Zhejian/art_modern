#include "art_profile_builder/exe/IntermediateEmpDist.hh"

#include "art_profile_builder/exe/IntermediateEmpDistPosition.hh"

#include "libam_support/utils/arithmetic_utils.hh"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace labw::art_modern {

IntermediateEmpDist::IntermediateEmpDist(const std::size_t read_length)
    : read_length_(read_length)
{
    positions_.resize(read_length);
}

bool IntermediateEmpDist::parse_read(const bam1_t* b)
{
    const auto* seq = bam_get_seq(b);
    const auto* qual = bam_get_qual(b);
    if (qual[0] == 0xff) {
        // No quality
        return false;
    }
    for (std::size_t i = 0; i < am_min(static_cast<std::size_t>(b->core.l_qseq), read_length_); ++i) {
        positions_[i].add(seq_nt16_str[bam_seqi(seq, i)], static_cast<am_qual_t>(qual[i]));
    }
    return true;
}

void IntermediateEmpDist::accumulate()
{
    for (auto& pos : positions_) {
        pos.accumulate();
    }
}

void IntermediateEmpDist::add(const IntermediateEmpDist& other)
{
    if (read_length_ != other.read_length_) {
        throw std::runtime_error("Incompatible IntermediateEmpDist sizes");
    }
    for (std::size_t i = 0; i < read_length_; ++i) {
        positions_[i].add(other.positions_[i]);
    }
}

void IntermediateEmpDist::write(std::ostream& oss) const
{
    for (const auto& base_idx : IntermediateEmpDistPosition::BASE_IDX) {
        for (std::size_t pos_id = 0; pos_id < positions_.size(); ++pos_id) {
            positions_[pos_id].write(oss, pos_id, base_idx);
        }
    }
}
} // namespace labw::art_modern
