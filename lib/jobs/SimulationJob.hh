#pragma once

#include "fasta/BaseFastaFetch.hh"
#include "fasta/CoverageInfo.hh"
#include "fasta/InMemoryFastaFetch.hh"

namespace labw {
namespace art_modern {

    class SimulationJob {
    public:
        ~SimulationJob();
        SimulationJob(BaseFastaFetch* fasta_fetch, const CoverageInfo& coverage_info, const int job_id);

        SimulationJob(SimulationJob&& other) noexcept:
        fasta_fetch(other.fasta_fetch),
        coverage_info(other.coverage_info),
        job_id(other.job_id)
        {
        };
        SimulationJob(const SimulationJob & ) = delete;
        SimulationJob& operator=(SimulationJob&& ) =delete;

        BaseFastaFetch* fasta_fetch;
        const CoverageInfo& coverage_info;
        const int job_id;
    };

} // art_modern
} // labw
