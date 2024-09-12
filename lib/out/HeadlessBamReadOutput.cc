#include "HeadlessBamReadOutput.hh"
#include "CExceptionsProxy.hh"
#include "global_variables.hh"
#include "seq_utils.hh"
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>

namespace po = boost::program_options;

namespace labw {
namespace art_modern {
    HeadlessBamReadOutput::HeadlessBamReadOutput(const std::string& filename, const SamOptions& sam_options)
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
        CExceptionsProxy::requires_numeric(
            sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
        BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
    }
    void HeadlessBamReadOutput::writeSE(const PairwiseAlignment& pwa)
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        if (is_closed_) {
            return;
        }
        auto sam_record = (bam1_t*)CExceptionsProxy::requires_not_null(
            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");
        auto rlen = static_cast<long>(pwa.query.size());
        auto seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);
        auto qual = pwa.qual;
        if (!pwa.is_plus_strand) {
            std::reverse(qual.begin(), qual.end());
        }
        auto oa_tag = bam_utils_.generate_oa_tag(pwa);
        auto tag_len = 2 + 1 + oa_tag.size() + 1;
        CExceptionsProxy::requires_numeric(bam_set1(sam_record, pwa.read_name.size(), pwa.read_name.c_str(),
                                               BAM_FUNMAP, // Alignment info moved to OA tag
                                               TID_FOR_UNMAPPED, // Alignment info moved to OA tag
                                               0, // Alignment info moved to OA tag
                                               0, // Alignment info moved to OA tag
                                               0, // Alignment info moved to OA tag
                                               nullptr, // Alignment info moved to OA tag
                                               TID_FOR_UNMAPPED, // Unset for SE reads
                                               0, // Unset for SE reads
                                               0, // Unset for SE reads
                                               rlen, seq.c_str(), qual.c_str(), tag_len),
            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        bam_aux_update_str(sam_record, "OA", static_cast<int>(oa_tag.length()), oa_tag.c_str());
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    }
    void HeadlessBamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        if (is_closed_) {
            return;
        }
        auto sam_record1 = (bam1_t*)CExceptionsProxy::requires_not_null(
            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");
        auto sam_record2 = (bam1_t*)CExceptionsProxy::requires_not_null(
            bam_init1(), USED_HTSLIB_NAME, "Failed to initialize SAM/BAM record");
        auto rlen = static_cast<long>(pwa1.query.size());
        auto seq1 = pwa1.is_plus_strand ? pwa1.query : revcomp(pwa1.query);
        auto seq2 = pwa2.is_plus_strand ? pwa2.query : revcomp(pwa2.query);
        auto qual1 = pwa1.qual;
        if (!pwa1.is_plus_strand) {
            std::reverse(qual1.begin(), qual1.end());
        }
        auto qual2 = pwa2.qual;
        if (!pwa1.is_plus_strand) {
            std::reverse(qual2.begin(), qual2.end());
        }
        auto oa_tag1 = bam_utils_.generate_oa_tag(pwa1);
        auto oa_tag2 = bam_utils_.generate_oa_tag(pwa2);
        auto tag_len1 = 2 + 1 + oa_tag1.size() + 1;
        auto tag_len2 = 2 + 1 + oa_tag2.size() + 1;
        CExceptionsProxy::requires_numeric(
            bam_set1(sam_record1, pwa1.read_name.size(), pwa1.read_name.c_str(),
                BAM_FPAIRED | BAM_FUNMAP | BAM_FMUNMAP | BAM_FREAD1, // Alignment info moved to OA tag
                TID_FOR_UNMAPPED, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                nullptr, // Alignment info moved to OA tag
                TID_FOR_UNMAPPED, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                rlen, seq1.c_str(), qual1.c_str(), tag_len1),
            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        bam_aux_update_str(sam_record1, "OA", static_cast<int>(oa_tag1.length()), oa_tag1.c_str());
        CExceptionsProxy::requires_numeric(
            bam_set1(sam_record2, pwa2.read_name.size(), pwa2.read_name.c_str(),
                BAM_FPAIRED | BAM_FUNMAP | BAM_FMUNMAP | BAM_FREAD2, // Alignment info moved to OA tag
                TID_FOR_UNMAPPED, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                nullptr, // Alignment info moved to OA tag
                TID_FOR_UNMAPPED, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                0, // Alignment info moved to OA tag
                rlen, seq2.c_str(), qual2.c_str(), tag_len2),
            USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        bam_aux_update_str(sam_record2, "OA", static_cast<int>(oa_tag2.length()), oa_tag2.c_str());

        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record1), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
        CExceptionsProxy::requires_numeric(sam_write1(sam_file_, sam_header_, sam_record2), USED_HTSLIB_NAME,
            "Failed to write SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    }
    void HeadlessBamReadOutput::close()
    {
        std::unique_lock<std::mutex> rhs_lk(mutex_);
        if (is_closed_) {
            return;
        }
        sam_close(sam_file_);
        is_closed_ = true;
    }
    HeadlessBamReadOutput::~HeadlessBamReadOutput() { HeadlessBamReadOutput::close(); }

    void HeadlessBamReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
    {
        po::options_description bam_desc("Headless SAM/BAM Output");
        bam_desc.add_options()("o-hl_sam", po::value<std::string>(),
            "Destination of output headless SAM/BAM file. Unset to disable the writer.");
        bam_desc.add_options()("o-hl_sam-use_m", "Whether to use CIGAR 'M' instead of '=/X' for alignment");
        bam_desc.add_options()("o-hl_sam-write_bam", "Enforce BAM instead of SAM output.");
        desc.add(bam_desc);
    }

    BaseReadOutput* HeadlessBamReadOutputFactory::create(
        const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch) const
    {
        if (vm.count("o-hl_sam")) {
            if (fasta_fetch->num_seqs() != 0) {
                BOOST_LOG_TRIVIAL(warning) << "Sequences presented in the reference file. Use SAM/BAM instead of this "
                                              "headless one for better compatibility.";
            }
            auto so = SamOptions();
            so.use_m = vm.count("o-hl_sam-use_m") > 0;
            so.write_bam = vm.count("o-hl_sam-write_bam") > 0;
            so.PG_CL = boost::algorithm::join(args, " ");
            return new HeadlessBamReadOutput(vm["o-hl_sam"].as<std::string>(), so);
        }
        return new DumbReadOutput();
    }

    HeadlessBamReadOutputFactory::~HeadlessBamReadOutputFactory() = default;
}
}
