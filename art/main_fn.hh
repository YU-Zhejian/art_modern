#pragma once
#include <exception>

#include "ArtContig.hh"
#include "ArtParams.hh"
#include "Empdist.hh"
#include "out/OutputDispatcher.hh"

namespace labw {
namespace art_modern {

    // when max_num =-1, no limit on the number of indels
    // the maxium number of indels is set by cdf_cutoff to save computation time
    void print_banner();

    void generate_all(const std::string& contig_name, const std::string& ref_seq, const ArtParams& art_params,
        const Empdist& qdist, double x_fold, OutputDispatcher& output_dispatcher);

    class GeneratedSeq {
    public:
        /**
         * The first strand FASTQ
         */
        std::string fastq;
        /**
         * The second strand FASTQ
         */
        std::string fastq2;
        /**
         * The strand SAM
         */
        std::string sam;
    };

} // namespace art_modern

} // namespace labw