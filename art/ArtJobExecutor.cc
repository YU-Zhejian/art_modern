#include "ArtJobExecutor.hh"
#include "ArtContig.hh"
#include "global_variables.hh"
#include "utils/mpi_utils.hh"
#include <boost/log/trivial.hpp>
#include <boost/timer/progress_display.hpp>

namespace labw::art_modern {

bool is_good(const ArtRead& art_read, const ArtParams& art_params, const std::string& qual_str)
{
    try {
        art_read.assess_num_n();
    } catch (TooMuchNException&) {
        return false;
    }
    if (art_read.seq_read.size() != art_params.read_len) {
        goto error;
    }
    if (qual_str.size() != art_params.read_len) {
        goto error;
    }
    if (art_read.aln_read.size() != art_read.aln_ref.size()) {
        goto error;
    }
    return true;
error:
    // #ifdef CEU_CM_IS_DEBUG
    //         abort_mpi();
    // #else
    return false; // FIXME: No idea why this occurs.
    // #endif
}

bool ArtJobExecutor::generate_pe(ArtContig& art_contig, const bool is_plus_strand)
{
    std::ostringstream osID;
    auto const& qdist = art_params_.qdist;

    std::vector<int> qual_1;
    std::vector<int> qual_2;

    osID << art_contig.seq_name_ << ':' << art_params_.id << ":" << read_id++;
    std::string sam_read_id = osID.str();
    ArtRead read_1(art_params_, rprob_);
    ArtRead read_2(art_params_, rprob_);

    try {
        art_contig.generate_read_pe(
            is_plus_strand, art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::MP, read_1, read_2);
    } catch (ReadGenerationException&) {
        return false;
    }

    if (!art_params_.sep_flag) {
        qdist.get_read_qual(qual_1, art_params_.read_len, rprob_, true);
        qdist.get_read_qual(qual_2, art_params_.read_len, rprob_, true);
    } else {
        qdist.get_read_qual_sep_1(qual_1, read_1.seq_read, rprob_);
        qdist.get_read_qual_sep_2(qual_2, read_1.seq_read, rprob_);
    }
    read_1.generate_snv_on_qual(qual_1);
    read_2.generate_snv_on_qual(qual_2);
    auto qual_1_str = qual_to_str(qual_1);
    auto qual_2_str = qual_to_str(qual_2);
    read_1.generate_pairwise_aln();
    read_2.generate_pairwise_aln();
    if (!(is_good(read_1, art_params_, qual_1_str) && is_good(read_2, art_params_, qual_2_str))) {
        return false;
    }

    output_dispatcher_->writePE(
        PairwiseAlignment(sam_read_id, art_contig.seq_name_, read_1.seq_read, read_1.seq_ref, qual_1_str,
            read_1.aln_read, read_1.aln_ref,
            read_1.is_plus_strand ? read_1.bpos : art_contig.ref_len_ - (read_1.bpos + art_params_.read_len),
            read_1.is_plus_strand),
        PairwiseAlignment(sam_read_id, art_contig.seq_name_, read_2.seq_read, read_2.seq_ref, qual_2_str,
            read_2.aln_read, read_2.aln_ref,
            read_2.is_plus_strand ? read_2.bpos : art_contig.ref_len_ - (read_2.bpos + art_params_.read_len),
            read_2.is_plus_strand));

    return true;
}

bool ArtJobExecutor::generate_se(ArtContig& art_contig, const bool is_plus_strand)

{
    std::ostringstream osID;
    std::vector<int> qual;
    std::string read_name;

    osID << art_contig.seq_name_ << ':' << art_params_.id << ':' << read_id++;
    read_name = osID.str();
    ArtRead art_read(art_params_, rprob_);
    try {
        art_contig.generate_read_se(is_plus_strand, art_read);
    } catch (ReadGenerationException&) {
        return false;
    }

    if (!art_params_.sep_flag) {
        art_params_.qdist.get_read_qual(qual, art_params_.read_len, rprob_, true);
    } else {
        art_params_.qdist.get_read_qual_sep_1(qual, art_read.seq_read, rprob_);
    }
    art_read.generate_snv_on_qual(qual);
    auto qual_str = qual_to_str(qual);

    art_read.generate_pairwise_aln();

    if (!is_good(art_read, art_params_, qual_str)) {
        return false;
    }

    output_dispatcher_->writeSE(PairwiseAlignment(read_name, art_contig.seq_name_, art_read.seq_read, art_read.seq_ref,
        qual_str, art_read.aln_read, art_read.aln_ref,
        art_read.is_plus_strand ? art_read.bpos : art_contig.ref_len_ - (art_read.bpos + art_params_.read_len),
        art_read.is_plus_strand));
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
    if (job_.fasta_fetch->num_seqs() == 0) {
        return;
    }
    const int num_reads_to_reduce = art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;

    BOOST_LOG_TRIVIAL(info) << "Starting simulation for job " << job_.job_id;
    for (size_t seq_id = 0; seq_id < job_.fasta_fetch->num_seqs(); ++seq_id) {
        const auto contig_name = job_.fasta_fetch->seq_name(seq_id);
        const auto contig_size = job_.fasta_fetch->seq_len(seq_id);
        if (job_.fasta_fetch->seq_len(seq_id) < art_params_.read_len) {
            BOOST_LOG_TRIVIAL(debug) << "Warning: the reference sequence " << contig_name << " (length "
                                     << job_.fasta_fetch->seq_len(seq_id)
                                     << "bps ) is skipped as it < the defined read length (" << art_params_.read_len
                                     << " bps)";
            continue;
        }
        BOOST_LOG_TRIVIAL(debug) << "Starting simulation for job " << job_.job_id << " CONTIG: " << contig_name;

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
            }
        }
        while (num_neg_reads > 0) {
            if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, false)
                                                                         : generate_se(art_contig, false)) {
                num_neg_reads -= num_reads_to_reduce;
            }
        }
    }
    BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id;
}

ArtJobExecutor::~ArtJobExecutor() = default;
}