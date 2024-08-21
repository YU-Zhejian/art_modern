#pragma once

#include "ArtParams.hh"
#include "Empdist.hh"
#include "out/OutputDispatcher.hh"

namespace labw {
namespace art_modern {
    // when max_num =-1, no limit on the number of indels
    // the maxium number of indels is set by cdf_cutoff to save computation time
    void print_banner();
    void handle_dumps();
    void init_logger();
    void generate_all(const ArtParams& art_params);

} // namespace art_modern

} // namespace labw