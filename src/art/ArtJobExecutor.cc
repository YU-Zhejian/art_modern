#include "ArtJobExecutor.hh"
#include "ArtConstants.hh"
#include "ArtContig.hh"
#include "utils/mpi_utils.hh"
#include <atomic>
#include <boost/log/trivial.hpp>
#include <boost/timer/progress_display.hpp>

namespace labw::art_modern {

void ArtJobExecutor::generate(const long targeted_num_reads, bool is_positive, ArtContig& art_contig)
{
    int num_cont_fail = 0;
    const auto max_tolerance
        = std::max(5L, static_cast<long>(static_cast<double>(targeted_num_reads) * MAX_TRIAL_RATIO_BEFORE_FAIL));
    const int num_reads_to_reduce = art_params.art_lib_const_mode == ART_LIB_CONST_MODE::SE ? 1 : 2;
    bool retv;
    long remaining_num_reads = targeted_num_reads;
    long current_num_reads = 0;
    while (remaining_num_reads > 0) {
        if (art_params.art_lib_const_mode == ART_LIB_CONST_MODE::SE) {
            retv = generate_se(art_contig, is_positive, current_num_reads);
        } else {
            retv = generate_pe(art_contig, is_positive, current_num_reads);
        }
        if (retv) {
            remaining_num_reads -= num_reads_to_reduce;
            num_reads += num_reads_to_reduce;
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

bool ArtJobExecutor::generate_pe(ArtContig& art_contig, const bool is_plus_strand, const std::size_t current_num_reads)
{
    std::ostringstream osID;
    osID << art_contig.seq_name << ':' << art_params.id << ":" << job_.job_id << ":" << mpi_rank_ << ":"
         << current_num_reads;
    const std::string read_name = osID.str();

    // std::snprintf(nullptr, 20, "%s:%s:%d:%s:%zu", art_params.id.c_str(), art_contig.seq_name.c_str(), job_.job_id,
    // mpi_rank_.c_str(),current_num_reads);
    ArtRead read_1(art_params, rprob_, art_contig.seq_name, read_name);
    ArtRead read_2(art_params, rprob_, art_contig.seq_name, read_name);

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

bool ArtJobExecutor::generate_se(ArtContig& art_contig, const bool is_plus_strand, const std::size_t current_num_reads)
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

        const auto num_pos_reads = static_cast<long>(cov_pos * cov_ratio);
        const auto num_neg_reads = static_cast<long>(cov_neg * cov_ratio);

        if (num_pos_reads + num_neg_reads == 0) {
            BOOST_LOG_TRIVIAL(debug) << "No read will be generated for the reference sequence " << contig_name
                                     << " due to insufficient coverage (pos=" << cov_pos << ", neg=" << cov_neg << ")";
            continue;
        }
        generate(num_pos_reads, true, art_contig);
        generate(num_neg_reads, false, art_contig);
    }

    BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id << " with " << num_reads
                            << " reads generated.";
    if (job_.free_fasta_fetch_after_execution) {
        delete job_.fasta_fetch;
    }
    is_running = false;
}
ArtJobExecutor::ArtJobExecutor(ArtJobExecutor&& other) noexcept
    : num_reads(other.num_reads.load())
    , art_params(other.art_params)
    , job_(std::move(other.job_))
    , rprob_(Rprob(art_params.pe_frag_dist_mean, art_params.pe_frag_dist_std_dev, art_params.read_len))
    , output_dispatcher_(other.output_dispatcher_)
    , mpi_rank_(other.mpi_rank_)

{
}
std::string ArtJobExecutor::thread_info() const { return std::to_string(job_.job_id) + ":" + mpi_rank_; }

ArtJobExecutor::~ArtJobExecutor() = default;
}