#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>

#include "BamReadOutput.hh"
#include "CExceptionsProxy.hh"
#include "global_variables.hh"
#include "seq_utils.hh"
namespace po = boost::program_options;

namespace labw {
namespace art_modern {

    void BamReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        if (is_closed_) {
            return;
        }

        int tid = CExceptionsProxy::requires_numeric(sam_hdr_name2tid(sam_header_, pwa.contig_name.c_str()),
            USED_HTSLIB_NAME, "Failed to fetch TID for contig '" + pwa.contig_name + "'", false,
            CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        auto sam_record = (bam1_t*)CExceptionsProxy::requires_not_null(
            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");
        auto rlen = static_cast<long>(pwa.query.size());
        auto cigar = pwa.generate_cigar_array(sam_options_.use_m);

        auto seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);
        auto qual = pwa.qual;
        hts_pos_t pos = pwa.align_contig_start;
        if (!pwa.is_plus_strand) {
            std::reverse(qual.begin(), qual.end());
        }
        auto cigar_c_arr = cigar_arr_to_c(cigar);
        assert_correct_cigar(pwa, cigar, cigar_c_arr);
        CExceptionsProxy::requires_numeric(
            bam_set1(sam_record, pwa.read_name.length(), pwa.read_name.c_str(), pwa.is_plus_strand ? 0 : BAM_FREVERSE,
                tid, pos, MAPQ_MAX, cigar.size(), cigar_c_arr,
                0, // Unset for SE reads
                0, // Unset for SE reads
                0, // Unset for SE reads
                rlen, seq.c_str(), qual.c_str(), 0),
            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        free(cigar_c_arr);
    }

    void BamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        if (is_closed_) {
            return;
        }

        int tid = CExceptionsProxy::requires_numeric(sam_hdr_name2tid(sam_header_, pwa1.contig_name.c_str()),
            USED_HTSLIB_NAME, "Failed to fetch TID for contig '" + pwa1.contig_name + "'", false,
            CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        auto rlen = static_cast<long>(pwa1.query.size());

        auto cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
        auto cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

        auto cigar1_arr = cigar_arr_to_c(cigar1);
        auto cigar2_arr = cigar_arr_to_c(cigar2);

        auto seq1 = pwa1.is_plus_strand ? pwa1.query : revcomp(pwa1.query);
        auto seq2 = pwa2.is_plus_strand ? pwa2.query : revcomp(pwa2.query);

        auto qual1 = pwa1.qual;
        auto qual2 = pwa2.qual;

        uint16_t flag1 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
        uint16_t flag2 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;

        hts_pos_t pos1;
        hts_pos_t pos2;

        pos1 = pwa1.align_contig_start;
        pos2 = pwa2.align_contig_start;
        if (pwa1.is_plus_strand) {
            flag1 |= BAM_FMREVERSE;
            flag2 |= BAM_FREVERSE;
            std::reverse(qual2.begin(), qual2.end());
        } else {
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

        auto sam_record1 = (bam1_t*)CExceptionsProxy::requires_not_null(
            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");
        auto sam_record2 = (bam1_t*)CExceptionsProxy::requires_not_null(
            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");

        assert_correct_cigar(pwa1, cigar1, cigar1_arr);
        assert_correct_cigar(pwa2, cigar2, cigar2_arr);

        CExceptionsProxy::requires_numeric(
            bam_set1(sam_record1, pwa1.read_name.length(), pwa1.read_name.c_str(), flag1, tid, pos1, MAPQ_MAX,
                cigar1.size(), cigar1_arr, tid, pos2, isize1, rlen, seq1.c_str(), qual1.c_str(), 0),
            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(
            bam_set1(sam_record2, pwa2.read_name.length(), pwa2.read_name.c_str(), flag2, tid, pos2, MAPQ_MAX,
                cigar2.size(), cigar2_arr, tid, pos1, isize2, rlen, seq2.c_str(), qual2.c_str(), 0),
            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record1), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record2), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

        free(cigar1_arr);
        free(cigar2_arr);
    }
    BamReadOutput::~BamReadOutput() { BamReadOutput::close(); }
    BamReadOutput::BamReadOutput(
        const std::string& filename, BaseFastaFetch* fasta_fetch, const SamReadOutputOptions& sam_options)
        : sam_options_(sam_options)
        , bam_utils_(sam_options)
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);

        sam_file_ = (samFile*)CExceptionsProxy::requires_not_null(
            sam_open(filename.c_str(), sam_options_.write_bam ? "wb" : "wh"), USED_HTSLIB_NAME,
            "Failed to open SAM file");
        sam_header_ = (sam_hdr_t*)CExceptionsProxy::requires_not_null(
            sam_hdr_init(), USED_HTSLIB_NAME, "Faield to initialize SAM header");
        CExceptionsProxy::requires_numeric(sam_hdr_add_line(sam_header_, "HD", "VN", sam_options_.HD_VN.c_str(), "SO",
                                               sam_options_.HD_SO.c_str(), NULL),
            USED_HTSLIB_NAME, "Failed to add HD header line");
        CExceptionsProxy::requires_numeric(sam_hdr_add_line(sam_header_, "PG", "ID", sam_options_.PG_ID.c_str(), "PN",
                                               sam_options_.PG_PN.c_str(), "CL", sam_options_.PG_CL.c_str(), NULL),
            USED_HTSLIB_NAME, "Failed to add PG header line", false, CExceptionsProxy::EXPECTATION::ZERO);

        fasta_fetch->update_sam_header(sam_header_);
        CExceptionsProxy::requires_numeric(
            sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
        BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
    }

    void BamReadOutput::close()
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        if (is_closed_) {
            return;
        }
        sam_close(sam_file_);
        is_closed_ = true;
    }
    void BamReadOutputFactory::patch_options(boost::program_options::options_description& desc)
    {
        po::options_description bam_desc("SAM/BAM Output");
        bam_desc.add_options()(
            "o-sam", po::value<std::string>(), "Destination of output SAM/BAM file. Unset to disable the writer.");
        bam_desc.add_options()("o-sam-use_m", "Whether to use CIGAR 'M' instead of '=/X' for alignment");
        bam_desc.add_options()("o-sam-write_bam", "Enforce BAM instead of SAM output.");
        desc.add(bam_desc);
    }
    BaseReadOutput* BamReadOutputFactory::create(
        const boost::program_options::variables_map& vm, BaseFastaFetch* fasta_fetch) const
    {
        if (vm.count("o-sam")) {
            if (fasta_fetch->num_seqs() == 0) {
                BOOST_LOG_TRIVIAL(error) << "No sequences in the reference file. If you used "
                                         << INPUT_FILE_PARSER_STREAM
                                         << " input parser, you should use headless SAM/BAM instead of this one.";
                exit(EXIT_FAILURE);
            }
            auto so = SamReadOutputOptions();
            so.use_m = vm.count("o-sam-use_m") > 0;
            so.write_bam = vm.count("o-sam-write_bam") > 0;
            so.PG_CL = boost::algorithm::join(args, " ");
            return new BamReadOutput(vm["o-sam"].as<std::string>(), fasta_fetch, so);
        }
        return new DumbReadOutput();
    }

    BamReadOutputFactory::~BamReadOutputFactory() = default;
}
}
