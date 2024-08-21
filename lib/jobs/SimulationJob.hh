#pragma once

#include "fasta/BaseFastaFetch.hh"
#include "fasta/CoverageInfo.hh"
#include <memory>

namespace labw {
namespace art_modern {

    class SimulationJob {
    public:
        SimulationJob(
            const std::shared_ptr<BaseFastaFetch>& fasta_fetch, const CoverageInfo& coverage_info, long job_id);
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch() const;
        const CoverageInfo& coverage_info() const;
        const long job_id;

    private:
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch_;
        const CoverageInfo& coverage_info_;
    };

} // art_modern
} // labw
