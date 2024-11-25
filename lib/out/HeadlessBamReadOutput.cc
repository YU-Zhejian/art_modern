#include "HeadlessBamReadOutput.hh"
#include "CExceptionsProxy.hh"
#include "art_modern_config.h"
#include "global_variables.hh"
#include "utils/seq_utils.hh"
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>

namespace po = boost::program_options;

namespace labw::art_modern {
HeadlessBamReadOutput::HeadlessBamReadOutput(const std::string& filename, const SamOptions& sam_options)
    : sam_options_(sam_options)
    , bam_utils_(sam_options)
{
    std::unique_lock rhs_lk(mutex_);
    sam_file_ = CExceptionsProxy::assert_not_null(
        sam_open(filename.c_str(), sam_options_.write_bam ? "wb" : "wh"), USED_HTSLIB_NAME, "Failed to open SAM file");
    sam_header_
        = CExceptionsProxy::assert_not_null(sam_hdr_init(), USED_HTSLIB_NAME, "Faield to initialize SAM header");
    CExceptionsProxy::assert_numeric(
        sam_hdr_add_line(sam_header_, "HD", "VN", sam_options_.HD_VN.c_str(), "SO", sam_options_.HD_SO.c_str(), NULL),
        USED_HTSLIB_NAME, "Failed to add HD header line");
    CExceptionsProxy::assert_numeric(sam_hdr_add_line(sam_header_, "PG", "ID", sam_options_.PG_ID.c_str(), "PN",
                                         sam_options_.PG_PN.c_str(), "CL", sam_options_.PG_CL.c_str(), NULL),
        USED_HTSLIB_NAME, "Failed to add PG header line", false, CExceptionsProxy::EXPECTATION::ZERO);
    CExceptionsProxy::assert_numeric(
        sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
}
void HeadlessBamReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    auto sam_record = BamUtils::init();
    const auto rlen = static_cast<long>(pwa.query.size());
    const auto& seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);
    auto qual = pwa.qual;
    auto cigar = pwa.generate_cigar_array(sam_options_.use_m);
    assert_correct_cigar(pwa, cigar);
    if (!pwa.is_plus_strand) {
        std::reverse(qual.begin(), qual.end());
        std::reverse(cigar.begin(), cigar.end());
    }

    const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa, cigar);
    const auto& oa_tag = BamUtils::generate_oa_tag(pwa, cigar, nm_tag);
    BamTags tags;
    tags.add_string("OA", oa_tag);
    tags.add_string("MD", md_tag);
    tags.add_int_i("NM", nm_tag);

    CExceptionsProxy::assert_numeric(bam_set1(sam_record, pwa.read_name.size(), pwa.read_name.c_str(),
                                         BAM_FUNMAP, // Alignment info moved to OA tag
                                         TID_FOR_UNMAPPED, // Alignment info moved to OA tag
                                         0, // Alignment info moved to OA tag
                                         0, // Alignment info moved to OA tag
                                         0, // Alignment info moved to OA tag
                                         nullptr, // Alignment info moved to OA tag
                                         TID_FOR_UNMAPPED, // Unset for SE reads
                                         0, // Unset for SE reads
                                         0, // Unset for SE reads
                                         rlen, seq.c_str(), qual.c_str(), tags.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    tags.patch(sam_record);

    std::unique_lock rhs_lk(mutex_);
    BamUtils::write(sam_file_, sam_header_, sam_record);
    bam_destroy1(sam_record);
}
void HeadlessBamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    auto sam_record1 = BamUtils::init();
    auto sam_record2 = BamUtils::init();

    const auto rlen = static_cast<long>(pwa1.query.size());

    const auto& seq1 = pwa1.is_plus_strand ? pwa1.query : revcomp(pwa1.query);
    const auto& seq2 = pwa2.is_plus_strand ? pwa2.query : revcomp(pwa2.query);

    auto qual1 = pwa1.qual;
    auto qual2 = pwa2.qual;

    auto cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
    auto cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

    assert_correct_cigar(pwa1, cigar1);
    assert_correct_cigar(pwa2, cigar2);

    if (!pwa1.is_plus_strand) {
        std::reverse(qual1.begin(), qual1.end());
        std::reverse(cigar1.begin(), cigar1.end());
    }
    if (!pwa1.is_plus_strand) {
        std::reverse(qual2.begin(), qual2.end());
        std::reverse(cigar2.begin(), cigar2.end());
    }

    const auto& [nm_tag1, md_tag1] = BamUtils::generate_nm_md_tag(pwa1, cigar1);
    const auto& [nm_tag2, md_tag2] = BamUtils::generate_nm_md_tag(pwa2, cigar2);

    const auto& oa_tag1 = BamUtils::generate_oa_tag(pwa1, cigar1, nm_tag1);
    const auto& oa_tag2 = BamUtils::generate_oa_tag(pwa2, cigar2, nm_tag2);

    BamTags tags1;
    tags1.add_string("OA", oa_tag1);
    tags1.add_string("MD", md_tag1);
    tags1.add_int_i("NM", nm_tag1);

    BamTags tags2;
    tags2.add_string("OA", oa_tag2);
    tags2.add_string("MD", md_tag2);
    tags2.add_int_i("NM", nm_tag2);

    CExceptionsProxy::assert_numeric(
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
            rlen, seq1.c_str(), qual1.c_str(), tags1.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(
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
            rlen, seq2.c_str(), qual2.c_str(), tags2.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    tags1.patch(sam_record1);
    tags2.patch(sam_record2);

    std::unique_lock rhs_lk(mutex_);
    BamUtils::write(sam_file_, sam_header_, sam_record1);
    BamUtils::write(sam_file_, sam_header_, sam_record2);
    bam_destroy1(sam_record1);
    bam_destroy1(sam_record2);
}
void HeadlessBamReadOutput::close()
{
    if (is_closed_) {
        return;
    }
    std::unique_lock rhs_lk(mutex_);
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

BaseReadOutput* HeadlessBamReadOutputFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const
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
