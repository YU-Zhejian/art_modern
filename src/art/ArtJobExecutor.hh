#pragma once
#include "art/ArtContig.hh"
#include "art/ArtParams.hh"
#include "art/random_generator.hh"

#include "libam/Dtypes.hh"
#include "libam/jobs/SimulationJob.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/utils/class_macros_utils.hh"

#include <atomic>
#include <cstddef>
#include <string>

namespace labw::art_modern {

class ArtJobExecutor {
public:
    ArtJobExecutor(ArtJobExecutor&& other) noexcept;

    DELETE_COPY(ArtJobExecutor)
    ArtJobExecutor& operator=(ArtJobExecutor&&) = delete;

    ArtJobExecutor(SimulationJob job, const ArtParams& art_params, BaseReadOutput* output_dispatcher);

    ~ArtJobExecutor() = default;
    void operator()();
    std::string thread_info() const;
    std::atomic<bool> is_running = false;

private:
    bool generate_pe(ArtContig& art_contig, bool is_plus_strand, std::size_t current_num_reads);
    bool generate_se(ArtContig& art_contig, bool is_plus_strand, std::size_t current_num_reads);
    void generate(am_readnum_t targeted_num_reads, bool is_positive, ArtContig& art_contig);

    const ArtParams& art_params_;
    SimulationJob job_;
    const std::string mpi_rank_;
    std::atomic<std::size_t> num_reads;
    BaseReadOutput* output_dispatcher_;
    Rprob rprob_;
    const int num_reads_to_reduce_;
};

} // namespace labw::art_modern
