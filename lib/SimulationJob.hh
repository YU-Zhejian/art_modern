#pragma once

#include "fasta/BaseFastaFetch.hh"
#include <memory>

namespace labw {
namespace art_modern {

    enum class SIMULATION_FRAGMENTATION_TYPE { TEMPLATE, RANDOM };

    class SimulationJob {
    public:
        SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, int num_reads_positive,
            int num_reads_negative, SIMULATION_FRAGMENTATION_TYPE fragmentation_type);
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch() const;
        int num_reads_positive() const;
        int num_reads_negative() const;
        SIMULATION_FRAGMENTATION_TYPE fragmentation_type() const;

    private:
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch_;
        int num_reads_positive_;
        int num_reads_negative_;
        SIMULATION_FRAGMENTATION_TYPE fragmentation_type_;
    };

} // art_modern
} // labw
