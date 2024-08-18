#pragma once

#include "fasta/BaseFastaFetch.hh"
#include "fasta/CoverageInfo.hh"
#include <memory>

namespace labw {
namespace art_modern {

    enum class SIMULATION_FRAGMENTATION_TYPE { TEMPLATE, RANDOM };

    class SimulationJob {
    public:
        SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, CoverageInfo coverage_info,
            SIMULATION_FRAGMENTATION_TYPE fragmentation_type);
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch() const;
        SIMULATION_FRAGMENTATION_TYPE fragmentation_type() const;
        const CoverageInfo& coverage_info() const;

    private:
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch_;
        CoverageInfo coverage_info_;
        SIMULATION_FRAGMENTATION_TYPE fragmentation_type_;
    };

} // art_modern
} // labw
