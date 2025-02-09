#pragma once

#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <memory>

namespace labw::art_modern {

class SimulationJob {
public:
    SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
        const std::shared_ptr<CoverageInfo>& coverage_info, int job_id);
    SimulationJob& operator=(SimulationJob&&) = delete;

    SimulationJob(SimulationJob&& other) noexcept;
    DELETE_COPY(SimulationJob)
    ~SimulationJob() = default;

    std::shared_ptr<BaseFastaFetch> fasta_fetch;
    std::shared_ptr<CoverageInfo> coverage_info;
    const int job_id;
};

} // namespace labw::art_modern
