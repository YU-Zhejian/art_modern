#include "ArtJobExecutor.hh"
#include "ArtConstants.hh"
#include "ArtContig.hh"
#include "utils/mpi_utils.hh"
#include <atomic>
#include <boost/log/trivial.hpp>
#include <boost/timer/progress_display.hpp>
#include <chrono>
#include <thread>

#define GENERATE(NUM_WHAT_READS, IS_POSITIVE, GEN_FUNC)                                                                \
    max_tolerance = std::max(5, static_cast<int>(NUM_WHAT_READS * MAX_TRIAL_RATIO_BEFORE_FAIL));                       \
    while (NUM_WHAT_READS > 0) {                                                                                       \
        if (GEN_FUNC(art_contig, IS_POSITIVE, num_reads)) {                                                            \
            NUM_WHAT_READS -= num_reads_to_reduce;                                                                     \
            num_reads += num_reads_to_reduce;                                                                          \
        } else {                                                                                                       \
            num_cont_fail++;                                                                                           \
        }                                                                                                              \
        if (num_cont_fail >= max_tolerance) {                                                                          \
            BOOST_LOG_TRIVIAL(debug) << "Failed to generate reads for " << contig_name << " sized " << contig_size     \
                                     << " after " << max_tolerance << " attempts.";                                    \
            break;                                                                                                     \
        }                                                                                                              \
    }

namespace labw::art_modern {

bool ArtJobExecutor::generate_pe(ArtContig& art_contig, bool is_plus_strand, const std::size_t current_num_reads)
{
    std::ostringstream osID;
    osID << art_contig.seq_name << ':' << art_params.id << ":" << job_.job_id << ":" << mpi_rank_ << ":"
         << current_num_reads;

    // std::snprintf(nullptr, 20, "%s:%s:%d:%s:%zu", art_params.id.c_str(), art_contig.seq_name.c_str(), job_.job_id, mpi_rank_.c_str(),current_num_reads);
    ArtRead read_1(art_params, rprob_, art_contig.seq_name, osID.str());
    ArtRead read_2(art_params, rprob_, art_contig.seq_name, osID.str());

    try {
        art_contig.generate_read_pe(
            is_plus_strand, art_params.art_lib_const_mode == ART_LIB_CONST_MODE::MP, read_1, read_2);
    } catch (ReadGenerationException&) {
        return false;
    }

    read_1.generate_snv_on_qual(true);
    read_2.generate_snv_on_qual(false);
    read_1.generate_pairwise_aln();
    read_2.generate_pairwise_aln();
    if (!(read_1.is_good() && read_2.is_good())) {
        return false;
    }

    output_dispatcher_->writePE(read_1.to_pwa(), read_2.to_pwa());

    return true;
}

bool ArtJobExecutor::generate_se(ArtContig& art_contig, bool is_plus_strand, const std::size_t current_num_reads)
{
    std::ostringstream osID;
    osID << art_contig.seq_name << ':' << art_params.id << ":" << job_.job_id << ":" << mpi_rank_ << ":"
         << current_num_reads;
    ArtRead art_read(art_params, rprob_, art_contig.seq_name, osID.str());
    try {
        art_contig.generate_read_se(is_plus_strand, art_read);
    } catch (ReadGenerationException&) {
        return false;
    }
    art_read.generate_snv_on_qual(true);
    art_read.generate_pairwise_aln();
    if (!art_read.is_good()) {
        return false;
    }
    output_dispatcher_->writeSE(art_read.to_pwa());
    return true;
}

ArtJobExecutor::ArtJobExecutor(SimulationJob job, const ArtParams& art_params, BaseReadOutput* output_dispatcher)
    : art_params(art_params)
    , job_(std::move(job))
    , rprob_(art_params.pe_frag_dist_mean, art_params.pe_frag_dist_std_dev, art_params.read_len)
    , output_dispatcher_(output_dispatcher)
    , mpi_rank_(mpi_rank())
{
    num_reads = 0;
}

void ArtJobExecutor::execute()
{
    is_running = true;
    const auto num_contigs = job_.fasta_fetch->num_seqs();
    if (num_contigs == 0) {
        return;
    }
    const int num_reads_to_reduce = art_params.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;

    BOOST_LOG_TRIVIAL(info) << "Starting simulation for job " << job_.job_id << " with " << num_contigs << " contigs";
    for (decltype(job_.fasta_fetch->num_seqs()) seq_id = 0; seq_id < num_contigs; ++seq_id) {
        const auto& contig_name = job_.fasta_fetch->seq_name(seq_id);
        const auto contig_size = job_.fasta_fetch->seq_len(seq_id);
        if (contig_size < art_params.read_len) {
            BOOST_LOG_TRIVIAL(debug) << "the reference sequence " << contig_name << " (length "
                                     << job_.fasta_fetch->seq_len(seq_id)
                                     << "bps ) is skipped as it < the defined read length (" << art_params.read_len
                                     << " bps)";
            continue;
        }

        ArtContig art_contig(job_.fasta_fetch, seq_id, art_params, rprob_);
        const double cov_ratio = art_params.art_simulation_mode == SIMULATION_MODE::TEMPLATE
            ? 1
            : static_cast<double>(contig_size) / art_params.read_len;
        const auto cov_pos = job_.coverage_info.coverage_positive(contig_name);
        const auto cov_neg = job_.coverage_info.coverage_negative(contig_name);

        auto num_pos_reads = static_cast<long>(cov_pos * cov_ratio);
        auto num_neg_reads = static_cast<long>(cov_neg * cov_ratio);

        if (num_pos_reads + num_neg_reads == 0) {
            BOOST_LOG_TRIVIAL(debug) << "No read will be generated for the reference sequence " << contig_name
                                     << " due to insufficient coverage (pos=" << cov_pos << ", neg=" << cov_neg << ")";
            continue;
        }
        int num_cont_fail = 0;
        int max_tolerance;

        if (art_params.art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
            GENERATE(num_pos_reads, true, generate_se)
            GENERATE(num_neg_reads, false, generate_se)
        } else {
            GENERATE(num_pos_reads, true, generate_pe)
            GENERATE(num_neg_reads, false, generate_pe)
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id << " with " << num_reads
                            << " reads generated.";
    if (job_.free_fasta_fetch_after_execution) {
        delete job_.fasta_fetch;
    }
    is_running = false;
}
ArtJobExecutor::ArtJobExecutor(ArtJobExecutor&& other) noexcept
    : art_params(other.art_params)
    , job_(std::move(other.job_))
    , rprob_(Rprob(art_params.pe_frag_dist_mean, art_params.pe_frag_dist_std_dev, art_params.read_len))
    , output_dispatcher_(other.output_dispatcher_)
    , mpi_rank_(other.mpi_rank_)
    , num_reads(other.num_reads.load())
{
}
std::string ArtJobExecutor::thread_info() const { return std::to_string(job_.job_id) + ":" + mpi_rank_; }

ArtJobExecutor::~ArtJobExecutor() = default;
}