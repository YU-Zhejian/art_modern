//
// Created by yuzj on 24-8-14.
//

#include "SimulationJob.hh"

namespace labw {
namespace art_modern {
    SimulationJob::SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch, const int num_reads_positive,
        const int num_reads_negative, const SIMULATION_FRAGMENTATION_TYPE fragmentation_type)
        : fasta_fetch_(fasta_fetch)
        , num_reads_positive_(num_reads_positive)
        , num_reads_negative_(num_reads_negative)
        , fragmentation_type_(fragmentation_type)
    {
    }
    const std::shared_ptr<BaseFastaFetch>& SimulationJob::fasta_fetch() const { return fasta_fetch_; }
    int SimulationJob::num_reads_positive() const { return num_reads_positive_; }
    int SimulationJob::num_reads_negative() const { return num_reads_negative_; }
    SIMULATION_FRAGMENTATION_TYPE SimulationJob::fragmentation_type() const { return fragmentation_type_; }
} // art_modern
} // labw