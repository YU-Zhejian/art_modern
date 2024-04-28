#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/asio.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>

#include "ArtContig.hh"
#include "ArtSamHeader.hh"
#include "SamRead.hh"
#include "ThreadSafeFileStream.hh"
#include "main_fn.hh"
#include "art_modern_constants.hh"
#include "seq_utils.hh"

using namespace std;

namespace labw {
namespace art_modern {

    void print_banner()
    {
        cout << endl;
        cout << "YuZJ Modified ART_Illumina" << endl;
        cout << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)" << endl;
        cout << "Originally written by: Weichun Huang <whduke@gmail.com>" << endl;
    }

    GeneratedSeq generate_pe(const ArtParams& art_params, const Empdist& qdist, const string& id, const ArtContig& art_contig, const long& t_num_read)
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

        auto arp = art_params.is_mp ? art_contig.generate_read_mp(art_params.is_amplicon) : art_contig.generate_read_pe(art_params.is_amplicon);
        if (art_params.mask_n) {
            if (arp.read_1.is_plus_strand) {
                size_t bpos2 = art_contig._ref_seq.size() - arp.read_2.bpos - art_params.read_len;
                if (art_contig._masked_pos.count(arp.read_1.bpos) > 0 || art_contig._masked_pos.count(bpos2) > 0) {
                    return generated_seq;
                }
            } else {
                size_t bpos1 = art_contig._ref_seq.size() - arp.read_1.bpos - art_params.read_len;
                if (art_contig._masked_pos.count(bpos1) > 0 || art_contig._masked_pos.count(arp.read_2.bpos) > 0) {
                    return generated_seq;
                }
            }
        }

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
            sam_read_1.tLen = static_cast<int>(sam_read_2.pos + arp.read_2.seq_read.size() - sam_read_1.pos);
            sam_read_2.tLen = -sam_read_1.tLen;
        } else {
            sam_read_2.tLen = static_cast<int>(sam_read_1.pos + arp.read_1.seq_read.size() - sam_read_2.pos);
            sam_read_1.tLen = -sam_read_2.tLen;
        }
        sam_read_1.printRead(SAMFILE);
        sam_read_2.printRead(SAMFILE);
        generated_seq.sam.append(SAMFILE.str());
        return generated_seq;
    }

    GeneratedSeq generate_se(const ArtParams& art_params, const Empdist& qdist, const string& id, const ArtContig& art_contig, const long& t_num_read)
    {
        ostringstream osID;
        ostringstream FQFILE;
        GeneratedSeq generated_seq;
        vector<int> qual;
        string read_id;

        osID << id << ':' << art_params.id << t_num_read;
        read_id = osID.str();
        auto art_read = art_contig.generate_read_se(art_params.is_amplicon);
        if (art_params.mask_n) {
            if (art_read.is_plus_strand) {
                if (art_contig._masked_pos.count(art_read.bpos) > 0) {
                    return generated_seq;
                }
            } else {
                size_t bpos = art_contig._ref_seq.size() - art_read.bpos - art_params.read_len;
                if (art_contig._masked_pos.count(bpos) > 0) {
                    return generated_seq;
                }
            }
        }

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
            sam_read.pos = art_contig._ref_seq.size() - (art_read.bpos + art_read.seq_read.size() - 1);
            sam_read.reverse_comp();
        }
        sam_read.printRead(SAMFILE);
        generated_seq.sam.append(SAMFILE.str());
        return generated_seq;
    }

    string generate_sam_header(const string& id, const string& ref_seq, const ArtParams& art_params)
    {
        vector<string> sn;
        vector<int> size;
        sn.emplace_back(id);
        size.emplace_back(ref_seq.size());
        ArtSamHeader sH(sn, size);

        sH.CL = boost::algorithm::join(art_params._args, " ");
        return sH.printHeader();
    }

    void generate_all(
        const string& id, const string& ref_seq, const ArtParams& art_params, const Empdist& qdist, double x_fold)
    {
        // TODO: no_sam not implemented!
        if (x_fold <= 0.0) {
            return;
        }
        ArtContig art_contig(ref_seq, id, art_params);
        if (art_contig._ref_seq.size() < art_params.read_len) {
            BOOST_LOG_TRIVIAL(warning) << "Warning: the reference sequence " << id << " (length " << art_contig._ref_seq.size() << "bps ) is skipped as it < the defined read length (" << art_params.read_len << " bps)" << endl;
            return;
        }
        ThreadSafeFileStream SAMFILE(art_params.samfile(id));
        if (!art_params.no_sam) {
            SAMFILE.write(generate_sam_header(id, ref_seq, art_params));
        }
        ThreadSafeFileStream FQFILE1(art_params.fqfile1(id));
        ThreadSafeFileStream FQFILE2(art_params.fqfile2(id));

        if (art_params.mask_n) {
            art_contig.mask_n_region(art_params.max_num_n);
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
        auto func = [num_read_per_batch, art_params, qdist, id, art_contig, &FQFILE1, &FQFILE2, &SAMFILE]() {
            auto current_num_read_per_batch = num_read_per_batch;
            while (current_num_read_per_batch > 0) {
                auto retv = (art_params.is_pe ? generate_pe : generate_se)(art_params, qdist, id, art_contig, current_num_read_per_batch);
                if (!retv.fastq.empty()) {
                    FQFILE1.write(retv.fastq);
                    FQFILE2.write(retv.fastq2);
                    if (!art_params.no_sam) {
                        SAMFILE.write(retv.sam);
                    }
                    current_num_read_per_batch -= art_params.is_pe ? 2 : 1;
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
        FQFILE2.close();
        SAMFILE.close();
    }
}
}
