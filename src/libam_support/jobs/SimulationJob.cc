#include "libam_support/jobs/SimulationJob.hh"

#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <memory>
#include <utility>

namespace labw::art_modern {

SimulationJob::SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
    const std::shared_ptr<CoverageInfo>& coverage_info, const int job_id)
    : fasta_fetch(fasta_fetch)
    , coverage_info(coverage_info)
    , job_id(job_id)
{
}
SimulationJob::SimulationJob(SimulationJob&& other) noexcept
    : fasta_fetch(std::move(other.fasta_fetch))
    , coverage_info(std::move(other.coverage_info))
    , job_id(other.job_id)
{
}
} // namespace labw::art_modern
