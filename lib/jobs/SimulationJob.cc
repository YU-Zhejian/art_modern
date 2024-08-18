#include <utility>

#include "SimulationJob.hh"

namespace labw {
namespace art_modern {
    const std::shared_ptr<BaseFastaFetch>& SimulationJob::fasta_fetch() const { return fasta_fetch_; }
    SIMULATION_FRAGMENTATION_TYPE SimulationJob::fragmentation_type() const { return fragmentation_type_; }
    SimulationJob::SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, CoverageInfo coverage_info,
        SIMULATION_FRAGMENTATION_TYPE fragmentation_type)
        : fasta_fetch_(fasta_fetch)
        , coverage_info_(std::move(coverage_info))
        , fragmentation_type_(fragmentation_type)
    {
    }
    const CoverageInfo& SimulationJob::coverage_info() const { return coverage_info_; }
} // art_modern
} // labw