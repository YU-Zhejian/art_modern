#pragma once
#include "ArtContig.hh"
#include "ArtParams.hh"
#include "jobs/SimulationJob.hh"
#include "out/BaseReadOutput.hh"
#include "random_generator.hh"

namespace labw {
namespace art_modern {

    class ArtJobExecutor {
    public:
        ArtJobExecutor(SimulationJob job, ArtParams art_params);
        void execute();

    private:
        bool generate_pe(ArtContig& art_contig, bool is_plus_strand);
        bool generate_se(ArtContig& art_contig, bool is_plus_strand);

        SimulationJob job_;
        ArtParams art_params_;
        Rprob rprob_;
        std::shared_ptr<BaseReadOutput> output_dispatcher_;
    };

} // art_modern
} // labw
