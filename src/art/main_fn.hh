#pragma once

#include "art/ArtParams.hh"

namespace labw::art_modern {
// when max_num =-1, no limit on the number of indels
// the maxium number of indels is set by cdf_cutoff to save computation time
void print_banner();
void generate_all(const ArtParams& art_params);

} // namespace labw::art_modern
