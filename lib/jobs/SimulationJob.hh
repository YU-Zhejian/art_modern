#pragma once

#include "fasta/BaseFastaFetch.hh"
#include "fasta/CoverageInfo.hh"
#include <memory>

namespace labw {
namespace art_modern {

    enum class SIMULATION_FRAGMENTATION_TYPE { TEMPLATE, RANDOM };

    class SimulationJob {
    public:
        SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, CoverageInfo coverage_info);
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch() const;
        const CoverageInfo& coverage_info() const;

    private:
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch_;
        CoverageInfo coverage_info_;
    };

} // art_modern
} // labw
