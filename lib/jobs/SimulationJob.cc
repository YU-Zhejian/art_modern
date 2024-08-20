#include <utility>

#include "SimulationJob.hh"

namespace labw {
namespace art_modern {
    const std::shared_ptr<BaseFastaFetch>& SimulationJob::fasta_fetch() const { return fasta_fetch_; }
    SimulationJob::SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, CoverageInfo coverage_info)
        : fasta_fetch_(fasta_fetch)
        , coverage_info_(std::move(coverage_info))
    {
    }
    const CoverageInfo& SimulationJob::coverage_info() const { return coverage_info_; }
} // art_modern
} // labw