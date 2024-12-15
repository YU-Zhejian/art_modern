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
    void operator()();
    std::string thread_info() const;
    std::atomic<bool> is_running = false;

private:
    bool generate_pe(ArtContig& art_contig, bool is_plus_strand, std::size_t current_num_reads);
    bool generate_se(ArtContig& art_contig, bool is_plus_strand, std::size_t current_num_reads);
    void generate(long targeted_num_reads, bool is_positive, ArtContig& art_contig);

    const ArtParams& art_params_;
    SimulationJob job_;
    const std::string mpi_rank_;
    std::atomic<std::size_t> num_reads;
    BaseReadOutput* output_dispatcher_;
    Rprob rprob_;
    std::vector<double> probs_indel_;
    const int num_reads_to_reduce;
};

} // namespace labw::art_modern
