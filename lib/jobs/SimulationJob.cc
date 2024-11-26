#include "SimulationJob.hh"

namespace labw::art_modern {
SimulationJob::~SimulationJob() =default;

SimulationJob::SimulationJob(BaseFastaFetch* fasta_fetch, const CoverageInfo& coverage_info, const int job_id, bool free_fasta_fetch_after_execution)
    : fasta_fetch(fasta_fetch)
    , coverage_info(coverage_info)
    , job_id(job_id)
, free_fasta_fetch_after_execution(free_fasta_fetch_after_execution)
{
}
SimulationJob::SimulationJob(SimulationJob&& other) noexcept
    : fasta_fetch(other.fasta_fetch)
    , coverage_info(other.coverage_info)
    , job_id(other.job_id)
    , free_fasta_fetch_after_execution(other.free_fasta_fetch_after_execution) {}
}