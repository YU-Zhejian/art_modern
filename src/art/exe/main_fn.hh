#pragma once

#include "art/lib/ArtIOParams.hh"
#include "art/lib/ArtParams.hh"

namespace labw::art_modern {
void print_banner();
void generate_all(const ArtParams& art_params, const ArtIOParams& art_io_params);

} // namespace labw::art_modern
