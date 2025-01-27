#include "art/ArtJobExecutor.hh"

#include "art/ArtConstants.hh"
#include "art/ArtContig.hh"
#include "art/ArtRead.hh"

#include "libam/Constants.hh"
#include "libam/Dtypes.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/utils/arithmetic_utils.hh"
#include "libam/utils/mpi_utils.hh"

#include <fmt/core.h>

#include <boost/log/trivial.hpp>

#include <htslib/hts.h>

#include <atomic>
#include <chrono>
#include <memory>
#include <string>
#include <thread>
#include <utility>

namespace labw::art_modern {
namespace {

    class AJEReporter {
    public:
        explicit AJEReporter(const ArtJobExecutor& aje)
            : aje_(aje)
        {
        }
        void stop()
        {
            should_stop_ = true;
            thread_.join();
        }
        void start() { thread_ = std::thread(&AJEReporter::job_, this); }

    private:
        void job_()
        {
            am_readnum_t num_reads_now = 0;
            am_readnum_t num_reads_prev = 0;
            std::this_thread::sleep_for(std::chrono::seconds(1));
            while (!should_stop_) {
                num_reads_now = aje_.num_reads();
                BOOST_LOG_TRIVIAL(info)
                    << "AJEReporter: Job " << aje_.thread_info() << ": " << num_reads_now << " reads generated. Speed: " << (num_reads_now - num_reads_prev) / 1 << "reads/s";
                num_reads_prev = num_reads_now;
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
        }

        const ArtJobExecutor& aje_;
        std::atomic<bool> should_stop_ { false };
        std::thread thread_;
    };

} // namespace
void ArtJobExecutor::generate(const am_readnum_t targeted_num_reads, const bool is_positive, ArtContig& art_contig)
{
    int num_cont_fail = 0;
    const auto max_tolerance = am_max(static_cast<am_readnum_t>(5),
        static_cast<am_readnum_t>(static_cast<double>(targeted_num_reads) * MAX_TRIAL_RATIO_BEFORE_FAIL));
    bool retv = false;
    am_readnum_t remaining_num_reads = targeted_num_reads;
    am_readnum_t current_num_reads = 0;
    while (remaining_num_reads > 0) {
        if (art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
            retv = generate_se(art_contig, is_positive, current_num_reads);
        } else {
            retv = generate_pe(art_contig, is_positive, current_num_reads);
        }
        if (retv) {
            remaining_num_reads -= num_reads_to_reduce_;
            num_reads_ += num_reads_to_reduce_;
            current_num_reads++;
        } else {
            num_cont_fail++;
        }
        if (num_cont_fail >= max_tolerance) {
            BOOST_LOG_TRIVIAL(debug) << "Failed to generate reads for " << art_contig.seq_name << " sized "
                                     << art_contig.seq_size << " after " << max_tolerance << " attempts.";
            break;
        }
    }
}

bool ArtJobExecutor::generate_pe(ArtContig& art_contig, const bool is_plus_strand, const am_readnum_t current_num_reads)
{
    const std::string read_name
        = fmt::format("{}:{}:{}:{}:{}", art_contig.seq_name, art_params_.id, job_.job_id, mpi_rank_, current_num_reads);

    ArtRead read_1(art_params_, art_contig.seq_name, read_name, rprob_);
    ArtRead read_2(art_params_, art_contig.seq_name, read_name, rprob_);
    art_contig.generate_read_pe(
        is_plus_strand, art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::MP, read_1, read_2);
    read_1.generate_snv_on_qual(true);
    read_2.generate_snv_on_qual(false);
    if (require_alignment_) {
        read_1.generate_pairwise_aln();
        read_2.generate_pairwise_aln();
    }
    if (!(read_1.is_good() && read_2.is_good())) {
        return false;
    }

    output_dispatcher_->writePE(read_1.to_pwa(), read_2.to_pwa());

    return true;
}

bool ArtJobExecutor::generate_se(ArtContig& art_contig, const bool is_plus_strand, const am_readnum_t current_num_reads)
{
    ArtRead art_read(art_params_, art_contig.seq_name,
        fmt::format("{}:{}:{}:{}:{}", art_contig.seq_name, art_params_.id, job_.job_id, mpi_rank_, current_num_reads),
        rprob_);
    art_contig.generate_read_se(is_plus_strand, art_read);
    art_read.generate_snv_on_qual(true);
    if (require_alignment_) {
        art_read.generate_pairwise_aln();
    }
    if (!art_read.is_good()) {
        return false;
    }
    output_dispatcher_->writeSE(art_read.to_pwa());
    return true;
}

ArtJobExecutor::ArtJobExecutor(
    SimulationJob&& job, const ArtParams& art_params, const std::shared_ptr<BaseReadOutput>& output_dispatcher)
    : art_params_(art_params)
    , job_(std::move(job))
    , mpi_rank_(mpi_rank())
    , num_reads_(0)
    , output_dispatcher_(output_dispatcher)
    , rprob_(art_params.pe_frag_dist_mean, art_params.pe_frag_dist_std_dev, art_params.read_len)
    , num_reads_to_reduce_(art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::SE ? 1 : 2)
    , require_alignment_(output_dispatcher->require_alignment())
{
}

bool ArtJobExecutor::is_running() const { return is_running_; }

am_readnum_t ArtJobExecutor::num_reads() const { return num_reads_; }

void ArtJobExecutor::operator()()
{
    is_running_ = true;
    AJEReporter reporter(*this);
    reporter.start();

    const auto num_contigs = job_.fasta_fetch->num_seqs();
    if (num_contigs == 0) {
        return;
    }
    hts_pos_t accumulated_contig_len = 0;

    BOOST_LOG_TRIVIAL(info) << "Starting simulation for job " << job_.job_id << " with " << num_contigs << " contigs";
    for (decltype(job_.fasta_fetch->num_seqs()) seq_id = 0; seq_id < num_contigs; ++seq_id) {
        const auto& contig_name = job_.fasta_fetch->seq_name(seq_id);
        const auto contig_size = job_.fasta_fetch->seq_len(seq_id);
        if (contig_size < art_params_.read_len) {
            BOOST_LOG_TRIVIAL(debug) << "the reference sequence " << contig_name << " (length "
                                     << job_.fasta_fetch->seq_len(seq_id)
                                     << "bps ) is skipped as it < the defined read length (" << art_params_.read_len
                                     << " bps)";
            continue;
        }
        accumulated_contig_len += contig_size;

        ArtContig art_contig(job_.fasta_fetch, seq_id, art_params_, rprob_);
        const double cov_ratio = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
            ? 1
            : static_cast<double>(contig_size) / art_params_.read_len;
        const auto cov_pos = job_.coverage_info->coverage_positive(contig_name);
        const auto cov_neg = job_.coverage_info->coverage_negative(contig_name);

        const auto num_pos_reads = static_cast<am_readnum_t>(cov_pos * cov_ratio);
        const auto num_neg_reads = static_cast<am_readnum_t>(cov_neg * cov_ratio);

        if (num_pos_reads + num_neg_reads == 0) {
            BOOST_LOG_TRIVIAL(debug) << "No read will be generated for the reference sequence " << contig_name
                                     << " due to insufficient coverage (pos=" << cov_pos << ", neg=" << cov_neg << ")";
            continue;
        }
        generate(num_pos_reads, true, art_contig);
        generate(num_neg_reads, false, art_contig);
    }

    BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id << " with " << num_reads_
                            << " reads (mean depth="
                            << static_cast<double>(num_reads_) * art_params_.read_len
            / static_cast<double>(accumulated_contig_len)
                            << ") generated.";
    reporter.stop();
    is_running_ = false;
}
ArtJobExecutor::ArtJobExecutor(ArtJobExecutor&& other) noexcept
    : art_params_(other.art_params_)
    , job_(std::move(other.job_))
    , mpi_rank_(other.mpi_rank_)
    , num_reads_(other.num_reads_.load())
    , output_dispatcher_(std::move(other.output_dispatcher_))
    , rprob_(Rprob(art_params_.pe_frag_dist_mean, art_params_.pe_frag_dist_std_dev, art_params_.read_len))
    , num_reads_to_reduce_(art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::SE ? 1 : 2)
    , require_alignment_(other.require_alignment_)
{
}
std::string ArtJobExecutor::thread_info() const { return std::to_string(job_.job_id) + ":" + mpi_rank_; }

} // namespace labw::art_modern
