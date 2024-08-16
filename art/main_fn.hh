#pragma once

#include "ArtParams.hh"
#include "Empdist.hh"
#include "out/OutputDispatcher.hh"

namespace labw {
namespace art_modern {
    // when max_num =-1, no limit on the number of indels
    // the maxium number of indels is set by cdf_cutoff to save computation time
    void print_banner();
    void generate_all(const std::string& contig_name, const std::string& ref_seq, const ArtParams& art_params,
        const Empdist& qdist, double x_fold, const std::shared_ptr<BaseReadOutput>& output_dispatcher);

} // namespace art_modern

} // namespace labw