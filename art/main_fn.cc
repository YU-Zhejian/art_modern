#include <atomic>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/asio.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>

#include <htslib/sam.h>

#include "ArtContig.hh"
#include "ArtSamHeader.hh"
#include "SamRead.hh"
#include "art_modern_constants.hh"
#include "main_fn.hh"
#include "seq_utils.hh"
#include "stream/DumbFileStream.hh"
#include "stream/ThreadSafeFileStream.hh"

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

    GeneratedSeq generate_pe(const ArtParams& art_params, const Empdist& qdist,
        const string& id, const ArtContig& art_contig,
        const long& t_num_read)
    {
        ostringstream osID;
        GeneratedSeq generated_seq;

        ostringstream FQFILE_1;
        ostringstream FQFILE_2;
        vector<int> qual_1;
        vector<int> qual_2;
        string read_id_1;
        string read_id_2;

        osID << id << ':' << art_params.id << t_num_read;
        string sam_read_id = osID.str();
        read_id_1 = sam_read_id + "/1";
        read_id_2 = sam_read_id + "/2";

        auto arp = art_params.art_lib_const_mode == ART_LIB_CONST_MODE::MP
            ? art_contig.generate_read_mp()
            : art_contig.generate_read_pe();

        if (!art_params.sep_flag) {
            qual_1 = qdist.get_read_qual(art_params.read_len, true);
            qual_2 = qdist.get_read_qual(art_params.read_len, false);
        } else {
            qual_1 = qdist.get_read_qual_sep_1(arp.read_1.seq_read);
            qual_2 = qdist.get_read_qual_sep_2(arp.read_2.seq_read);
        }

        qual_1 = arp.read_1.generate_snv_on_qual(qual_1);
        qual_2 = arp.read_2.generate_snv_on_qual(qual_2);

        FQFILE_1 << "@" << read_id_1 << endl
                 << arp.read_1.seq_read << endl
                 << "+" << endl
                 << qual_to_str(qual_1) << endl;
        generated_seq.fastq.append(FQFILE_1.str());

        auto aln_1 = arp.read_1.generate_pairwise_aln();

        FQFILE_2 << "@" << read_id_2 << endl
                 << arp.read_2.seq_read << endl
                 << "+" << endl
                 << qual_to_str(qual_2) << endl;
        auto aln_2 = arp.read_2.generate_pairwise_aln();
        generated_seq.fastq2.append(FQFILE_2.str());
        if (art_params.no_sam) {
            return generated_seq;
        }

        ostringstream SAMFILE;
        SamRead sam_read_1;
        SamRead sam_read_2;
        sam_read_1.rNext = "=";
        sam_read_2.rNext = "=";

        sam_read_1.qname = sam_read_id;
        sam_read_1.rname = id;
        sam_read_1.seq = arp.read_1.seq_read;
        sam_read_1.qual = qual_to_str(qual_1);

        sam_read_2.qname = sam_read_id;
        sam_read_2.rname = id;
        sam_read_2.seq = arp.read_2.seq_read;
        sam_read_2.qual = qual_to_str(qual_2);

        sam_read_1.flag = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
        sam_read_2.flag = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;
        if (arp.read_1.is_plus_strand) {
            sam_read_1.pos = arp.read_1.bpos + 1;
            sam_read_1.flag = sam_read_1.flag | BAM_FMREVERSE;
            sam_read_2.flag = sam_read_2.flag | BAM_FREVERSE;
            sam_read_2.pos = art_contig._ref_seq.size() - (arp.read_2.bpos + art_params.read_len - 1);
            sam_read_2.reverse_comp();
        } else {
            sam_read_1.pos = art_contig._ref_seq.size() - (arp.read_1.bpos + art_params.read_len - 1);
            sam_read_1.reverse_comp();
            sam_read_2.pos = arp.read_2.bpos + 1;
            sam_read_1.flag = sam_read_1.flag | BAM_FREVERSE;
            sam_read_2.flag = sam_read_2.flag | BAM_FMREVERSE;
        }
        sam_read_1.cigar = aln_1.generate_cigar(!arp.read_1.is_plus_strand, art_params.cigar_use_m);
        sam_read_2.cigar = aln_2.generate_cigar(!arp.read_2.is_plus_strand, art_params.cigar_use_m);

        sam_read_1.pNext = sam_read_2.pos;
        sam_read_2.pNext = sam_read_1.pos;
        if (sam_read_2.pos > sam_read_1.pos) {
            sam_read_1.tLen = static_cast<int>(
                sam_read_2.pos + art_params.read_len - sam_read_1.pos);
            sam_read_2.tLen = -sam_read_1.tLen;
        } else {
            sam_read_2.tLen = static_cast<int>(
                sam_read_1.pos + art_params.read_len - sam_read_2.pos);
            sam_read_1.tLen = -sam_read_2.tLen;
        }
        sam_read_1.printRead(SAMFILE);
        sam_read_2.printRead(SAMFILE);
        generated_seq.sam.append(SAMFILE.str());
        return generated_seq;
    }

    GeneratedSeq generate_se(const ArtParams& art_params, const Empdist& qdist,
        const string& id, const ArtContig& art_contig,
        const long& t_num_read)
    {
        ostringstream osID;
        ostringstream FQFILE;
        GeneratedSeq generated_seq;
        vector<int> qual;
        string read_id;

        osID << id << ':' << art_params.id << t_num_read;
        read_id = osID.str();
        auto art_read = art_contig.generate_read_se();

        if (!art_params.sep_flag) {
            qual = qdist.get_read_qual(art_params.read_len, true);
        } else {
            qual = qdist.get_read_qual_sep_1(art_read.seq_read);
        }
        qual = art_read.generate_snv_on_qual(qual);

        FQFILE << "@" << read_id << endl
               << art_read.seq_read << endl
               << "+" << endl
               << qual_to_str(qual) << endl;
        generated_seq.fastq.append(FQFILE.str());
        if (art_params.no_sam) {
            return generated_seq;
        }

        SamRead sam_read;
        ostringstream SAMFILE;

        auto aln = art_read.generate_pairwise_aln();
        sam_read.cigar = aln.generate_cigar(!art_read.is_plus_strand, art_params.cigar_use_m);

        sam_read.qname = read_id;
        sam_read.rname = id;
        sam_read.flag = 0;
        sam_read.seq = art_read.seq_read;
        sam_read.qual = qual_to_str(qual);
        if (art_read.is_plus_strand) {
            sam_read.pos = art_read.bpos + 1;
        } else {
            sam_read.flag = BAM_FREVERSE;
            sam_read.pos = art_contig._ref_seq.size() - (art_read.bpos + art_params.read_len - 1);
            sam_read.reverse_comp();
        }
        sam_read.printRead(SAMFILE);
        generated_seq.sam.append(SAMFILE.str());
        return generated_seq;
    }

    string generate_sam_header(const string& id, const string& ref_seq,
        const ArtParams& art_params)
    {
        vector<string> sn;
        vector<int> size;
        sn.emplace_back(id);
        size.emplace_back(ref_seq.size());
        ArtSamHeader sH(sn, size);

        sH.CL = boost::algorithm::join(art_params._args, " ");
        return sH.printHeader();
    }

    void generate_all(const string& id, const string& ref_seq,
        const ArtParams& art_params, const Empdist& qdist,
        double x_fold)
    {
        if (x_fold <= 0.0) {
            return;
        }
        ArtContig art_contig(ref_seq, id, art_params);
        if (art_contig._ref_seq.size() < art_params.read_len) {
            BOOST_LOG_TRIVIAL(warning)
                << "Warning: the reference sequence " << id << " (length "
                << art_contig._ref_seq.size()
                << "bps ) is skipped as it < the defined read length ("
                << art_params.read_len << " bps)";
            return;
        }
        std::shared_ptr<FileStreamInterface> SAMFILE_ptr;
        if (art_params.no_sam) {
            SAMFILE_ptr = std::make_shared<DumbFileStream>();
        } else {
            SAMFILE_ptr = std::make_shared<ThreadSafeFileStream>(art_params.samfile(id));
        }

        SAMFILE_ptr->write(generate_sam_header(id, ref_seq, art_params));

        ThreadSafeFileStream FQFILE1(art_params.fqfile1(id));

        std::shared_ptr<FileStreamInterface> FQFILE2_ptr;
        if (art_params.art_lib_const_mode != ART_LIB_CONST_MODE::SE) {
            FQFILE2_ptr = std::make_shared<ThreadSafeFileStream>(art_params.fqfile2(id));
        } else {
            FQFILE2_ptr = std::make_shared<DumbFileStream>();
        }

        auto t_num_read = static_cast<long>(static_cast<double>(art_contig._ref_seq.size()) / art_params.read_len * x_fold);
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
        auto func = [num_read_per_batch, art_params, qdist, id, art_contig, &FQFILE1,
                        FQFILE2_ptr, SAMFILE_ptr, &read_id]() {
            auto current_num_read_per_batch = num_read_per_batch;
            while (current_num_read_per_batch > 0) {
                auto retv = (art_params.art_lib_const_mode != ART_LIB_CONST_MODE::SE ? generate_pe : generate_se)(
                    art_params, qdist, id, art_contig, read_id++);
                if (!retv.fastq.empty()) {
                    FQFILE1.write(retv.fastq);
                    FQFILE2_ptr->write(retv.fastq2);
                    SAMFILE_ptr->write(retv.sam);
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

        FQFILE1.close();
        FQFILE2_ptr->close();
        SAMFILE_ptr->close();
    }
} // namespace art_modern
} // namespace labw
