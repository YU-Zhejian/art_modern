#pragma once

#include "fasta/BaseFastaFetch.hh"
#include <memory>

namespace labw {
namespace art_modern {

    enum class SIMULATION_FRAGMENTATION_TYPE { TEMPLATE, RANDOM };

    class SimulationJob {

    public:
        std::shared_ptr<BaseFastaFetch> fasta_fetch;
        int num_reads;
        SIMULATION_FRAGMENTATION_TYPE fragmentation_type;
    };

} // art_modern
} // labw
