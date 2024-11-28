#pragma once
#include "ArtContig.hh"
#include "ArtParams.hh"
#include "jobs/SimulationJob.hh"
#include "out/BaseReadOutput.hh"
#include "random_generator.hh"
#include <atomic>

namespace labw::art_modern {

class ArtJobExecutor {
public:
    ArtJobExecutor(ArtJobExecutor&& other) noexcept;
    ArtJobExecutor(const ArtJobExecutor&) = delete;
    ArtJobExecutor& operator=(ArtJobExecutor&&) = delete;
    ArtJobExecutor(SimulationJob job, const ArtParams& art_params, BaseReadOutput* output_dispatcher);

    ~ArtJobExecutor();
    void execute();
    std::atomic<std::size_t> num_reads;
    std::string thread_info() const;
    const ArtParams& art_params;
    std::atomic<bool> is_running = false;

private:
    bool generate_pe(ArtContig& art_contig, bool is_plus_strand, std::size_t current_num_reads);
    bool generate_se(ArtContig& art_contig, bool is_plus_strand, std::size_t current_num_reads);

    SimulationJob job_;
    Rprob rprob_;
    BaseReadOutput* output_dispatcher_;
    const std::string mpi_rank_;
};

} // namespace labw::art_modern
