#pragma once
#include "art/ArtContig.hh"
#include "art/ArtParams.hh"
#include "art/random_generator.hh"

#include "libam/Dtypes.hh"
#include "libam/jobs/JobExecutor.hh"
#include "libam/jobs/SimulationJob.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/utils/class_macros_utils.hh"

#include <atomic>
#include <memory>
#include <string>

namespace labw::art_modern {

class ArtJobExecutor : public JobExecutor {
public:
    ArtJobExecutor(
        SimulationJob&& job, const ArtParams& art_params, const std::shared_ptr<BaseReadOutput>& output_dispatcher);
    ~ArtJobExecutor() override = default;

    ArtJobExecutor(ArtJobExecutor&& other) noexcept;
    ArtJobExecutor& operator=(ArtJobExecutor&&) = delete;
    DELETE_COPY(ArtJobExecutor)
    [[nodiscard]] bool is_running() const override;

    void operator()() override;
    [[nodiscard]] std::string thread_info() const override;

private:
    bool generate_pe(ArtContig& art_contig, bool is_plus_strand, am_readnum_t current_num_reads);
    bool generate_se(ArtContig& art_contig, bool is_plus_strand, am_readnum_t current_num_reads);
    void generate(am_readnum_t targeted_num_reads, bool is_positive, ArtContig& art_contig);

    const ArtParams& art_params_;
    SimulationJob job_;
    const std::string mpi_rank_;
    am_readnum_t total_num_reads_generated_ = 0;
    std::shared_ptr<BaseReadOutput> output_dispatcher_;
    Rprob rprob_;
    const int num_reads_to_reduce_;
    const bool require_alignment_;
    std::atomic<bool> is_running_ = false;
    std::string current_contig_;
    am_readnum_t current_n_reads_left_ = 0;
    am_readnum_t current_n_fails_ = 0;
    am_readnum_t current_max_tolerence_ = 0;
    am_readnum_t current_n_reads_generated_ = 0;
};

} // namespace labw::art_modern
