#pragma once

#include "libam/ds/CoverageInfo.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

namespace labw::art_modern {

class SimulationJob {
public:
    ~SimulationJob();
    SimulationJob(BaseFastaFetch* fasta_fetch, const CoverageInfo& coverage_info, int job_id,
        bool free_fasta_fetch_after_execution);

    SimulationJob(SimulationJob&& other) noexcept;

    SimulationJob(const SimulationJob&) = delete;
    SimulationJob& operator=(SimulationJob&&) = delete;
    SimulationJob& operator=(const SimulationJob&) = delete;

    BaseFastaFetch* fasta_fetch;
    const CoverageInfo coverage_info; // We own this
    const int job_id;
    const bool free_fasta_fetch_after_execution;
};

} // namespace labw::art_modern
