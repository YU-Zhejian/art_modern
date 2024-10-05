#include "ArtJobExecutor.hh"
#include "ArtContig.hh"
#include "global_variables.hh"
#include <boost/log/trivial.hpp>
#include <boost/timer/progress_display.hpp>

namespace labw {
namespace art_modern {

    bool ArtJobExecutor::generate_pe(ArtContig& art_contig, const bool is_plus_strand)
    {
        std::ostringstream osID;
        auto const& qdist = art_params_.qdist;

        std::vector<int> qual_1;
        std::vector<int> qual_2;

        osID << art_contig.id_ << ':' << art_params_.id << ":" << read_id++;
        std::string sam_read_id = osID.str();
        ArtReadPair arp = art_params_.art_lib_const_mode == ART_LIB_CONST_MODE::MP
            ? art_contig.generate_read_mp(is_plus_strand)
            : art_contig.generate_read_pe(is_plus_strand);
        try {
            arp.read_1.assess_num_n();
            arp.read_2.assess_num_n();
        } catch (TooMuchNException&) {
            return false;
        }

        if (!art_params_.sep_flag) {
            qdist.get_read_qual(qual_1, art_params_.read_len, rprob_, true);
            qdist.get_read_qual(qual_2, art_params_.read_len, rprob_, true);
        } else {
            qdist.get_read_qual_sep_1(qual_1, arp.read_1.seq_read, rprob_);
            qdist.get_read_qual_sep_2(qual_2, arp.read_1.seq_read, rprob_);
        }

        auto qual_1_str = qual_to_str(arp.read_1.generate_snv_on_qual(qual_1));
        auto qual_2_str = qual_to_str(arp.read_2.generate_snv_on_qual(qual_2));

        arp.read_1.generate_pairwise_aln();
        arp.read_2.generate_pairwise_aln();

        if (arp.read_1.seq_read.size() != art_params_.read_len) {
            return false; // FIXME: No idea why this occurs.
        }

        if (arp.read_2.seq_read.size() != art_params_.read_len) {
            return false; // FIXME: No idea why this occurs.
        }

        output_dispatcher_->writePE(
            PairwiseAlignment(sam_read_id, art_contig.id_, arp.read_1.seq_read, arp.read_1.seq_ref, qual_1_str,
                arp.read_1.aln_read, arp.read_1.aln_ref,
                arp.read_1.is_plus_strand ? arp.read_1.bpos
                                          : art_contig.ref_len_ - (arp.read_1.bpos + art_params_.read_len),
                arp.read_1.is_plus_strand),
            PairwiseAlignment(sam_read_id, art_contig.id_, arp.read_2.seq_read, arp.read_2.seq_ref, qual_2_str,
                arp.read_2.aln_read, arp.read_2.aln_ref,
                arp.read_2.is_plus_strand ? arp.read_2.bpos
                                          : art_contig.ref_len_ - (arp.read_2.bpos + art_params_.read_len),
                arp.read_2.is_plus_strand));

        return true;
    }

    bool ArtJobExecutor::generate_se(ArtContig& art_contig, const bool is_plus_strand)

    {
        std::ostringstream osID;
        std::ostringstream FQFILE;
        std::vector<int> qual;
        std::string read_name;
        auto const& qdist = art_params_.qdist;

        osID << art_contig.id_ << ':' << art_params_.id << ":" << read_id++;
        read_name = osID.str();
        auto art_read = art_contig.generate_read_se(is_plus_strand);
        try {
            art_read.assess_num_n();
        } catch (TooMuchNException&) {
            return false;
        }
        if (!art_params_.sep_flag) {
            qdist.get_read_qual(qual, art_params_.read_len, rprob_, true);
        } else {
            qdist.get_read_qual_sep_1(qual, art_read.seq_read, rprob_);
        }
        auto qual_str = qual_to_str(art_read.generate_snv_on_qual(qual));
        art_read.generate_pairwise_aln();

        if (art_read.seq_read.size() != art_params_.read_len) {
            return false; // FIXME: No idea why this occurs.
        }

        output_dispatcher_->writeSE(PairwiseAlignment(read_name, art_contig.id_, art_read.seq_read, art_read.seq_ref,
            qual_str, art_read.aln_read, art_read.aln_ref,
            art_read.is_plus_strand ? art_read.bpos : art_contig.ref_len_ - (art_read.bpos + art_params_.read_len),
            art_read.is_plus_strand));
        return true;
    }

    ArtJobExecutor::ArtJobExecutor(SimulationJob job, const ArtParams& art_params)
        : job_(std::move(job))
        , art_params_(art_params)
        , rprob_(static_cast<float>(art_params_.pe_frag_dist_mean),
              static_cast<float>(art_params_.pe_frag_dist_std_dev), art_params_.read_len)
        , output_dispatcher_(art_params_.out_dispatcher)
    {
    }

    void ArtJobExecutor::execute()
    {
        if (job_.fasta_fetch->num_seqs() == 0) {
            return;
        }
        BOOST_LOG_TRIVIAL(info) << "Starting simulation for job " << job_.job_id;
        for (size_t seq_id = 0; seq_id < job_.fasta_fetch->num_seqs(); ++seq_id) {
            const std::string contig_name = job_.fasta_fetch->seq_name(seq_id);
            BOOST_LOG_TRIVIAL(debug) << "Starting simulation for job " << job_.job_id << " CONTIG: " << contig_name;

            ArtContig art_contig(job_.fasta_fetch, seq_id, art_params_, rprob_);

            BOOST_LOG_TRIVIAL(debug) << "Starting simulation for job " << job_.job_id << " CONTIG: " << contig_name
                                     << ": ArtContig created";
            if (art_contig.ref_len_ < art_params_.read_len) {
                // BOOST_LOG_TRIVIAL(warning)
                //     << "Warning: the reference sequence " << contig_name << " (length " << art_contig.ref_len_
                //     << "bps ) is skipped as it < the defined read length (" << art_params_.read_len << " bps)";
                BOOST_LOG_TRIVIAL(debug) << "Starting simulation for job " << job_.job_id << " CONTIG: " << contig_name
                                         << ": SKIPPED";
                continue;
            }
            auto coverage_positive = job_.coverage_info.coverage_positive(contig_name);
            auto coverage_negative = job_.coverage_info.coverage_negative(contig_name);

            long num_pos_reads;
            long num_neg_reads;

            if (art_params_.art_simulation_mode == SIMULATION_MODE::TEMPLATE) {
                num_pos_reads = static_cast<long>(coverage_positive);
                num_neg_reads = static_cast<long>(coverage_negative);
            } else {
                num_pos_reads = static_cast<long>(
                    coverage_positive * static_cast<double>(art_contig.ref_len_) / art_params_.read_len);
                num_neg_reads = static_cast<long>(
                    coverage_negative * static_cast<double>(art_contig.ref_len_) / art_params_.read_len);
            }

            while (num_pos_reads > 0) {
                BOOST_LOG_TRIVIAL(debug) << "Simulation for job " << job_.job_id << " CONTIG: " << contig_name
                                         << ": POS: " << num_pos_reads << " remaining";
                if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, true)
                                                                             : generate_se(art_contig, true)) {
                    num_pos_reads -= art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;
                }
            }
            while (num_neg_reads > 0) {
                BOOST_LOG_TRIVIAL(debug) << "Simulation for job " << job_.job_id << " CONTIG: " << contig_name
                                         << ": NEG: " << num_neg_reads << " remaining";
                if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, false)
                                                                             : generate_se(art_contig, false)) {
                    num_neg_reads -= art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;
                }
            }
        }
        BOOST_LOG_TRIVIAL(info) << "Finished simulation for job " << job_.job_id;
    }

    ArtJobExecutor::~ArtJobExecutor() = default;
} // art_modern
} // labw