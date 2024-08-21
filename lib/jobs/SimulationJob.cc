#include "ExceptionUtils.hh"

#include "SimulationJob.hh"

namespace labw {
namespace art_modern {
    const std::shared_ptr<BaseFastaFetch>& SimulationJob::fasta_fetch() const { return fasta_fetch_; }
    SimulationJob::SimulationJob(
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch, const CoverageInfo& coverage_info, const long job_id)
        : job_id(job_id)
        , fasta_fetch_(fasta_fetch)
        , coverage_info_(coverage_info)
    {
        if (fasta_fetch == nullptr) {
            throw_with_trace(std::runtime_error("fasta_fetch is nullptr"));
        }
    }
    const CoverageInfo& SimulationJob::coverage_info() const { return coverage_info_; }
} // art_modern
} // labw