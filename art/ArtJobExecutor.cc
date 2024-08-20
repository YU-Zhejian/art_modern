#include "ArtJobExecutor.hh"
#include "ArtContig.hh"
#include "global_variables.hh"
#include <boost/log/trivial.hpp>

using namespace std;

namespace labw {
namespace art_modern {

    bool ArtJobExecutor::generate_pe(ArtContig& art_contig, bool is_plus_strand)
    {
        ostringstream osID;
        auto const& qdist = art_params_.qdist;

        vector<int> qual_1;
        vector<int> qual_2;

        osID << art_contig.id_ << ':' << art_params_.id << read_id;
        string sam_read_id = osID.str();
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
            qual_1 = qdist.get_read_qual(art_params_.read_len, rprob_, true);
            qual_2 = qdist.get_read_qual(art_params_.read_len, rprob_, false);
        } else {
            qual_1 = qdist.get_read_qual_sep_1(arp.read_1.seq_read, rprob_);
            qual_2 = qdist.get_read_qual_sep_2(arp.read_2.seq_read, rprob_);
        }

        auto qual_1_str = qual_to_str(arp.read_1.generate_snv_on_qual(qual_1));
        auto qual_2_str = qual_to_str(arp.read_2.generate_snv_on_qual(qual_2));

        arp.read_1.generate_pairwise_aln();
        arp.read_2.generate_pairwise_aln();

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

    bool ArtJobExecutor::generate_se(ArtContig& art_contig, bool is_plus_strand)

    {
        ostringstream osID;
        ostringstream FQFILE;
        vector<int> qual;
        string read_name;
        auto const& qdist = art_params_.qdist;

        osID << art_contig.id_ << ':' << art_params_.id << read_id;
        read_name = osID.str();
        auto art_read = art_contig.generate_read_se(is_plus_strand);
        try {
            art_read.assess_num_n();
        } catch (TooMuchNException&) {
            return false;
        }
        if (!art_params_.sep_flag) {
            qual = qdist.get_read_qual(art_params_.read_len, rprob_, true);
        } else {
            qual = qdist.get_read_qual_sep_1(art_read.seq_read, rprob_);
        }
        auto qual_str = qual_to_str(art_read.generate_snv_on_qual(qual));
        art_read.generate_pairwise_aln();
        output_dispatcher_->writeSE(PairwiseAlignment(read_name, art_params_.id, art_read.seq_read, art_read.seq_ref,
            qual_str, art_read.aln_read, art_read.aln_ref,
            art_read.is_plus_strand ? art_read.bpos : art_contig.ref_len_ - (art_read.bpos + art_params_.read_len),
            art_read.is_plus_strand));
        return true;
    }

    ArtJobExecutor::ArtJobExecutor(SimulationJob job, ArtParams art_params)
        : job_(std::move(job))
        , art_params_(std::move(art_params))
        , rprob_(
              static_cast<float>(art_params_.pe_frag_dist_mean), static_cast<float>(art_params_.pe_frag_dist_std_dev))
        , output_dispatcher_(std::move(art_params.out_dispatcher))
    {
    }

    void ArtJobExecutor::execute()
    {
        for (const auto& contig_name : job_.fasta_fetch()->seq_names()) {

            ArtContig art_contig(job_.fasta_fetch(), contig_name, art_params_, rprob_);
            if (art_contig.ref_len_ < art_params_.read_len) {
                BOOST_LOG_TRIVIAL(warning)
                    << "Warning: the reference sequence " << contig_name << " (length " << art_contig.ref_len_
                    << "bps ) is skipped as it < the defined read length (" << art_params_.read_len << " bps)";
                return;
            }
            auto coverage_positive = job_.coverage_info().coverage_positive(contig_name);
            auto coverage_negative = job_.coverage_info().coverage_negative(contig_name);

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
                if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, true)
                                                                             : generate_se(art_contig, true)) {
                    num_pos_reads -= art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;
                }
            }
            while (num_neg_reads > 0) {
                if (art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe(art_contig, false)
                                                                             : generate_se(art_contig, false)) {
                    num_neg_reads -= art_params_.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;
                }
            }
        }
    }
} // art_modern
} // labw