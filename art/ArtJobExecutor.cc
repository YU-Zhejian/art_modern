#include "ArtJobExecutor.hh"
#include "ArtContig.hh"
#include "global_variables.hh"
#include "utils/mpi_utils.hh"
#include <boost/log/trivial.hpp>
#include <boost/timer/progress_display.hpp>

namespace labw::art_modern {

bool ArtJobExecutor::generate_pe(ArtContig& art_contig, const bool is_plus_strand)
{
    std::ostringstream osID;
    osID << art_contig.seq_name << ':' << art_params_.id << ":" << read_id++;
    ArtRead read_1(art_params_, rprob_, art_contig.seq_name, osID.str());
    ArtRead read_2(art_params_, rprob_, art_contig.seq_name, osID.str());

    try {
        art_contig.generate_read_pe(
            is_plus_strand, art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::MP, read_1, read_2);
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

bool ArtJobExecutor::generate_se(ArtContig& art_contig, const bool is_plus_strand)

{
    std::ostringstream osID;
    osID << art_contig.seq_name << ':' << art_params_.id << ':' << read_id++;
    ArtRead art_read(art_params_, rprob_, art_contig.seq_name, osID.str());
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

ArtJobExecutor::ArtJobExecutor(SimulationJob job, const ArtParams& art_params)
    : job_(std::move(job))
    , art_params_(art_params)
    , rprob_(art_params_.pe_frag_dist_mean, art_params_.pe_frag_dist_std_dev, art_params_.read_len)
    , output_dispatcher_(art_params_.out_dispatcher)
{
}

void ArtJobExecutor::execute()
{
    size_t num_reads = 0;
    const auto num_contigs = job_.fasta_fetch->num_seqs();
    if (num_contigs == 0) {
        return;
    }
    const int num_reads_to_reduce = art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;

    BOOST_LOG_TRIVIAL(info) << "Starting simulation for job " << job_.job_id << " with " << num_contigs << " contigs";
    for (auto seq_id = 0; seq_id < num_contigs; ++seq_id) {
        const auto contig_name = job_.fasta_fetch->seq_name(seq_id);
        const auto contig_size = job_.fasta_fetch->seq_len(seq_id);
        if (contig_size < art_params_.read_len) {
            BOOST_LOG_TRIVIAL(debug) << "the reference sequence " << contig_name << " (length "
                                     << job_.fasta_fetch->seq_len(seq_id)
                                     << "bps ) is skipped as it < the defined read length (" << art_params_.read_len
                                     << " bps)";
            continue;
        }
        ArtContig art_contig(job_.fasta_fetch, seq_id, art_params_, rprob_);
        const double cov_ratio = art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE
            ? 1
            : static_cast<double>(contig_size) / art_params_.read_len;

        auto num_pos_reads = static_cast<long>(job_.coverage_info.coverage_positive(contig_name) * cov_ratio);
        auto num_neg_reads = static_cast<long>(job_.coverage_info.coverage_negative(contig_name) * cov_ratio);

        while (num_pos_reads > 0) {
            if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, true)
                                                                         : generate_se(art_contig, true)) {
                num_pos_reads -= num_reads_to_reduce;
                num_reads += num_reads_to_reduce;
            }
        }
        while (num_neg_reads > 0) {
            if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, false)
                                                                         : generate_se(art_contig, false)) {
                num_neg_reads -= num_reads_to_reduce;
                num_reads += num_reads_to_reduce;
            }
        }
    }
    BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id << " with " << num_reads
                            << " reads generated.";
}

ArtJobExecutor::~ArtJobExecutor() = default;
}