#include <atomic>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/asio.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>

#include "ArtContig.hh"
#include "art_modern_constants.hh"
#include "main_fn.hh"
#include "seq_utils.hh"

using namespace std;

namespace labw {
namespace art_modern {

    void print_banner()
    {
        BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina";
        BOOST_LOG_TRIVIAL(info) << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
        BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
        BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <Zhejian.23@intl.zju.edu.cn>";
    }

    bool generate_pe(const ArtParams& art_params, const Empdist& qdist, const string& contig_name,
        const ArtContig& art_contig, const long& t_num_read, const std::shared_ptr<BaseReadOutput>& output_dispatcher)
    {
        ostringstream osID;

        vector<int> qual_1;
        vector<int> qual_2;

        osID << contig_name << ':' << art_params.id << t_num_read;
        string sam_read_id = osID.str();

        auto arp = art_params.art_lib_const_mode == ART_LIB_CONST_MODE::MP ? art_contig.generate_read_mp()
                                                                           : art_contig.generate_read_pe();

        if (!art_params.sep_flag) {
            qual_1 = qdist.get_read_qual(art_params.read_len, true);
            qual_2 = qdist.get_read_qual(art_params.read_len, false);
        } else {
            qual_1 = qdist.get_read_qual_sep_1(arp.read_1.seq_read);
            qual_2 = qdist.get_read_qual_sep_2(arp.read_2.seq_read);
        }

        auto qual_1_str = qual_to_str(arp.read_1.generate_snv_on_qual(qual_1));
        auto qual_2_str = qual_to_str(arp.read_2.generate_snv_on_qual(qual_2));

        arp.read_1.generate_pairwise_aln();
        arp.read_2.generate_pairwise_aln();

        output_dispatcher->writePE(PairwiseAlignment(sam_read_id, contig_name, arp.read_1.seq_read, arp.read_1.seq_ref,
                                       qual_1_str, arp.read_1.aln_read, arp.read_1.aln_ref,
                                       static_cast<hts_pos_t>(arp.read_1.is_plus_strand ? arp.read_1.bpos
                                                                                        : art_contig._ref_seq.length()
                                                   - (arp.read_1.bpos + art_params.read_len)),
                                       arp.read_1.is_plus_strand),
            PairwiseAlignment(sam_read_id, contig_name, arp.read_2.seq_read, arp.read_2.seq_ref, qual_2_str,
                arp.read_2.aln_read, arp.read_2.aln_ref,
                static_cast<hts_pos_t>(arp.read_2.is_plus_strand
                        ? arp.read_2.bpos
                        : art_contig._ref_seq.length() - (arp.read_2.bpos + art_params.read_len)),
                arp.read_2.is_plus_strand));

        return true;
    }

    bool generate_se(const ArtParams& art_params, const Empdist& qdist, const string& contig_name,
        const ArtContig& art_contig, const long& t_num_read, const std::shared_ptr<BaseReadOutput>& output_dispatcher)
    {
        ostringstream osID;
        ostringstream FQFILE;
        vector<int> qual;
        string read_id;

        osID << contig_name << ':' << art_params.id << t_num_read;
        read_id = osID.str();
        auto art_read = art_contig.generate_read_se();

        if (!art_params.sep_flag) {
            qual = qdist.get_read_qual(art_params.read_len, true);
        } else {
            qual = qdist.get_read_qual_sep_1(art_read.seq_read);
        }
        auto qual_str = qual_to_str(art_read.generate_snv_on_qual(qual));
        art_read.generate_pairwise_aln();
        output_dispatcher->writeSE(PairwiseAlignment(read_id, contig_name, art_read.seq_read, art_read.seq_ref,
            qual_str, art_read.aln_read, art_read.aln_ref,
            static_cast<hts_pos_t>(art_read.is_plus_strand
                    ? art_read.bpos
                    : art_contig._ref_seq.length() - (art_read.bpos + art_params.read_len)),
            art_read.is_plus_strand));
        return true;
    }

    void generate_all(const string& contig_name, const string& ref_seq, const ArtParams& art_params,
        const Empdist& qdist, double x_fold, const std::shared_ptr<BaseReadOutput>& output_dispatcher)
    {
        if (x_fold <= 0.0) {
            return;
        }
        ArtContig art_contig(ref_seq, contig_name, art_params);
        if (art_contig._ref_seq.size() < art_params.read_len) {
            BOOST_LOG_TRIVIAL(warning) << "Warning: the reference sequence " << contig_name << " (length "
                                       << art_contig._ref_seq.size()
                                       << "bps ) is skipped as it < the defined read length (" << art_params.read_len
                                       << " bps)";
            return;
        }

        auto t_num_read
            = static_cast<long>(static_cast<double>(art_contig._ref_seq.size()) / art_params.read_len * x_fold);
        int num_cores;

        if (art_params.parallel_on_read && art_params.parallel != PARALLEL_DISABLE) {
            if (art_params.parallel == PARALLEL_ALL) {
                num_cores = static_cast<int>(boost::thread::hardware_concurrency());
            } else {
                num_cores = art_params.parallel;
            }
        } else {
            num_cores = 1;
        }
        auto num_read_per_batch = static_cast<int>(t_num_read / num_cores + 1);
        std::atomic_long read_id;
        auto func
            = [&num_read_per_batch, &art_params, &qdist, &contig_name, &art_contig, &output_dispatcher, &read_id]() {
                  auto current_num_read_per_batch = num_read_per_batch;
                  while (current_num_read_per_batch > 0) {
                      if ((art_params.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe : generate_se)(
                              art_params, qdist, contig_name, art_contig, read_id++, output_dispatcher)) {
                          current_num_read_per_batch -= art_params.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? 2 : 1;
                      }
                  }
              };
        boost::asio::thread_pool pool(num_cores);
        for (int i = 0; i < num_cores; i++) {
            if (art_params.parallel_on_read && art_params.parallel != PARALLEL_DISABLE) {
                post(pool, func);
            } else {
                func();
            }
        }
        pool.join();
    }
} // namespace art_modern
} // namespace labw
