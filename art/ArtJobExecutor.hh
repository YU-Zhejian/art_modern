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
        ArtJobExecutor(ArtJobExecutor&& other) noexcept:
            job_(std::move(other.job_)),
            art_params_(other.art_params_),
            rprob_(other.rprob_),
            output_dispatcher_(other.output_dispatcher_){
        };
        ArtJobExecutor(const ArtJobExecutor & ) =delete;
        ArtJobExecutor& operator=(ArtJobExecutor&& ) =delete;

        ArtJobExecutor(SimulationJob job, const ArtParams& art_params);
        ~ArtJobExecutor();
        void execute();

    private:
        bool generate_pe(ArtContig& art_contig, bool is_plus_strand);
        bool generate_se(ArtContig& art_contig, bool is_plus_strand);

        SimulationJob job_;
        const ArtParams& art_params_;
        Rprob rprob_;
        BaseReadOutput* output_dispatcher_;
    };

} // art_modern
} // labw