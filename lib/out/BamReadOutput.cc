#include <algorithm>
#include <fstream>
#include <iostream>

#include "BamReadOutput.hh"
#include "CExceptionsProxy.hh"
#include "seq_utils.hh"

namespace labw {
namespace art_modern {
    void BamReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        int tid = CExceptionsProxy::requires_numeric(sam_hdr_name2tid(sam_header_, pwa.contig_name.c_str()),
            "htslib", UNKNOWN_C_EXCEPTION, false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        auto contig_len = sam_hdr_tid2len(sam_header_, tid);
        auto sam_record = (bam1_t*)CExceptionsProxy::requires_not_null(bam_init1(), "htslib");
        auto rlen = static_cast<long>(pwa.query.size());
        auto cigar = pwa.generate_cigar_array(!pwa.is_plus_strand, sam_options_.use_m);

        auto seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);
        auto qual = pwa.qual;
        hts_pos_t pos;
        if (pwa.is_plus_strand) {
            pos = pwa.align_contig_start + 1L;
        } else {
            pos = contig_len - (pwa.align_contig_start + rlen) + 1L;
            std::reverse(qual.begin(), qual.end());
        }
        CExceptionsProxy::requires_numeric(bam_set1(
                                               sam_record,
                                               pwa.read_name.length(),
                                               pwa.read_name.c_str(),
                                               pwa.is_plus_strand ? 0 : BAM_FREVERSE,
                                               tid,
                                               pos,
                                               MAPQ_MAX,
                                               cigar.size(),
                                               cigar.data(),
                                               0, // Unset for SE reads
                                               0, // Unset for SE reads
                                               0, // Unset for SE reads
                                               rlen,
                                               seq.c_str(),
                                               qual.c_str(),
                                               0),
            "htslib", UNKNOWN_C_EXCEPTION, false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record), "htslib");
    }

    void BamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        int tid = CExceptionsProxy::requires_numeric(sam_hdr_name2tid(sam_header_, pwa1.contig_name.c_str()),
            "htslib", UNKNOWN_C_EXCEPTION, false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        auto contig_len = sam_hdr_tid2len(sam_header_, tid);
        auto rlen = static_cast<long>(pwa1.query.size());

        auto cigar1 = pwa1.generate_cigar_array(!pwa1.is_plus_strand, sam_options_.use_m);
        auto cigar2 = pwa2.generate_cigar_array(!pwa2.is_plus_strand, sam_options_.use_m);

        auto seq1 = pwa1.is_plus_strand ? pwa1.query : revcomp(pwa1.query);
        auto seq2 = pwa2.is_plus_strand ? pwa2.query : revcomp(pwa2.query);

        auto qual1 = pwa1.qual;
        auto qual2 = pwa2.qual;

        uint16_t flag1 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
        uint16_t flag2 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;

        hts_pos_t pos1;
        hts_pos_t pos2;

        if (pwa1.is_plus_strand) {
            pos1 = pwa1.align_contig_start + 1L;
            pos2 = contig_len - (pwa2.align_contig_start + rlen) + 1L;
            flag1 |= BAM_FMREVERSE;
            flag2 |= BAM_FREVERSE;
            std::reverse(qual2.begin(), qual2.end());
        } else {
            pos1 = contig_len - (pwa1.align_contig_start + rlen) + 1L;
            pos2 = pwa2.align_contig_start + 1L;
            flag1 |= BAM_FREVERSE;
            flag2 |= BAM_FMREVERSE;
            std::reverse(qual1.begin(), qual1.end());
        }

        hts_pos_t isize1;
        hts_pos_t isize2;

        if (pos2 > pos1) {
            isize1 = pos2 + rlen - pos1;
            isize2 = -isize1;
        } else {
            isize2 = pos1 + rlen - pos2;
            isize1 = -isize2;
        }

        auto sam_record1 = (bam1_t*)CExceptionsProxy::requires_not_null(bam_init1(), "htslib");
        auto sam_record2 = (bam1_t*)CExceptionsProxy::requires_not_null(bam_init1(), "htslib");

        CExceptionsProxy::requires_numeric(bam_set1(
                                               sam_record1,
                                               pwa1.read_name.length(),
                                               pwa1.read_name.c_str(),
                                               flag1,
                                               tid,
                                               pos1,
                                               MAPQ_MAX,
                                               cigar1.size(),
                                               cigar1.data(),
                                               tid,
                                               pos2,
                                               isize1,
                                               rlen,
                                               seq1.c_str(),
                                               qual1.c_str(),
                                               0),
            "htslib", UNKNOWN_C_EXCEPTION, false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(bam_set1(
                                               sam_record2,
                                               pwa2.read_name.length(),
                                               pwa2.read_name.c_str(),
                                               flag2,
                                               tid,
                                               pos2,
                                               MAPQ_MAX,
                                               cigar2.size(),
                                               cigar2.data(),
                                               tid,
                                               pos1,
                                               isize2,
                                               rlen,
                                               seq2.c_str(),
                                               qual2.c_str(),
                                               0),
            "htslib", UNKNOWN_C_EXCEPTION, false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record1), "htslib");
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record2), "htslib");
    }
    BamReadOutput::~BamReadOutput()
    {
        sam_close(sam_file_);
    }
    BamReadOutput::BamReadOutput(const std::string& filename,
        const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
        const SamReadOutputOptions& sam_options)
        : sam_options_(sam_options)
    {
        sam_file_ = (samFile*)CExceptionsProxy::requires_not_null(sam_open(filename.c_str(), "wb"), "htslib");
        sam_header_ = (sam_hdr_t*)CExceptionsProxy::requires_not_null(sam_hdr_init(), "htslib");
        CExceptionsProxy::requires_numeric(
            sam_hdr_add_line(sam_header_, "HD", "HD_VN", sam_options_.HD_VN.c_str(), "HD_SO", sam_options_.HD_SO.c_str(), NULL),
            "htslib");
        CExceptionsProxy::requires_numeric(sam_hdr_add_line(
                                               sam_header_,
                                               "PG",
                                               "ID", sam_options_.PG_ID.c_str(),
                                               "PN", sam_options_.PG_PN.c_str(),
                                               "CL", sam_options_.PG_CL.c_str(),
                                               NULL),
            "htslib", UNKNOWN_C_EXCEPTION, false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

        fasta_fetch->update_sam_header(sam_header_);
        CExceptionsProxy::requires_numeric(sam_hdr_write(sam_file_, sam_header_),
            "htslib");
    }
}
}
