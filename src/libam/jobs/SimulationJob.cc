#include "libam/jobs/SimulationJob.hh"

#include "libam/ds/CoverageInfo.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

#include <memory>
#include <utility>

namespace labw::art_modern {

SimulationJob::SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
    const std::shared_ptr<CoverageInfo>& coverage_info, const int job_id, const bool free_fasta_fetch_after_execution)
    : fasta_fetch(fasta_fetch)
    , coverage_info(coverage_info)
    , job_id(job_id)
    , free_fasta_fetch_after_execution(free_fasta_fetch_after_execution)
{
}
SimulationJob::SimulationJob(SimulationJob&& other) noexcept
    : fasta_fetch(std::move(other.fasta_fetch))
    , coverage_info(std::move(other.coverage_info))
    , job_id(other.job_id)
    , free_fasta_fetch_after_execution(other.free_fasta_fetch_after_execution)
{
}
} // namespace labw::art_modern
