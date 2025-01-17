#include "art_modern_config.h" // For USED_HTSLIB_NAME

#include "libam/out/HeadlessBamReadOutput.hh"

#include "libam/CExceptionsProxy.hh"
#include "libam/Constants.hh"
#include "libam/bam/BamOptions.hh"
#include "libam/bam/BamTags.hh"
#include "libam/bam/BamUtils.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseFileReadOutput.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/DumbReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/utils/mpi_utils.hh"
#include "libam/utils/seq_utils.hh"

#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <htslib/sam.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace po = boost::program_options;

namespace labw::art_modern {
HeadlessBamReadOutput::HeadlessBamReadOutput(const std::string& filename, const BamOptions& sam_options)
    : BaseFileReadOutput(filename)
    , sam_file_(BamUtils::open_file(filename, sam_options))
    , sam_header_(BamUtils::init_header(sam_options))
    , sam_options_(sam_options)
    , lfio_(sam_file_, sam_header_)
{
    CExceptionsProxy::assert_numeric(
        sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
    lfio_.start();
}
void HeadlessBamReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    auto sam_record = BamUtils::init_uptr();
    const auto rlen = pwa.query.size();
    auto seq = pwa.query;
    auto qual = pwa.qual;
    auto cigar = pwa.generate_cigar_array(sam_options_.use_m);
    BamUtils::assert_correct_cigar(pwa, cigar);
    if (!pwa.is_plus_strand) {
        std::reverse(qual.begin(), qual.end());
        std::reverse(cigar.begin(), cigar.end());
        revcomp_inplace(seq);
    }

    const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa, cigar);
    const auto& oa_tag = BamUtils::generate_oa_tag(pwa, cigar, nm_tag);
    BamTags tags;
    tags.add_string("OA", oa_tag);
    tags.add_string("MD", md_tag);
    tags.add_int_i("NM", nm_tag);

    CExceptionsProxy::assert_numeric(bam_set1(sam_record.get(), pwa.read_name.size(), pwa.read_name.c_str(),
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
    tags.patch(sam_record.get());
    lfio_.push(std::move(sam_record));
}
void HeadlessBamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    auto sam_record1 = BamUtils::init_uptr();
    auto sam_record2 = BamUtils::init_uptr();

    const auto rlen = pwa1.query.size();

    auto seq1 = pwa1.query;
    auto seq2 = pwa2.query;

    auto cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
    auto cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

    BamUtils::assert_correct_cigar(pwa1, cigar1);
    BamUtils::assert_correct_cigar(pwa2, cigar2);

    if (!pwa1.is_plus_strand) {
        std::reverse(cigar1.begin(), cigar1.end());
        revcomp_inplace(seq1);
    }
    if (!pwa2.is_plus_strand) {
        std::reverse(cigar2.begin(), cigar2.end());
        revcomp_inplace(seq2);
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
        bam_set1(sam_record1.get(), pwa1.read_name.size(), pwa1.read_name.c_str(),
            BAM_FPAIRED | BAM_FUNMAP | BAM_FMUNMAP | BAM_FREAD1, // Alignment info moved to OA tag
            TID_FOR_UNMAPPED, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            nullptr, // Alignment info moved to OA tag
            TID_FOR_UNMAPPED, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            rlen, seq1.c_str(), pwa1.qual.c_str(), tags1.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record2.get(), pwa2.read_name.size(), pwa2.read_name.c_str(),
            BAM_FPAIRED | BAM_FUNMAP | BAM_FMUNMAP | BAM_FREAD2, // Alignment info moved to OA tag
            TID_FOR_UNMAPPED, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            nullptr, // Alignment info moved to OA tag
            TID_FOR_UNMAPPED, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            rlen, seq2.c_str(), pwa2.qual.c_str(), tags2.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    if (!pwa1.is_plus_strand) {
        reverse(bam_get_qual(sam_record1), rlen);
    } else {
        reverse(bam_get_qual(sam_record2), rlen);
    }

    tags1.patch(sam_record1.get());
    tags2.patch(sam_record2.get());

    lfio_.push(std::move(sam_record1));
    lfio_.push(std::move(sam_record2));
}
void HeadlessBamReadOutput::close()
{
    if (is_closed_) {
        return;
    }
    lfio_.stop();
    sam_close(sam_file_);
    sam_hdr_destroy(sam_header_);
    BaseFileReadOutput::close();
}
HeadlessBamReadOutput::~HeadlessBamReadOutput() { HeadlessBamReadOutput::close(); }

void HeadlessBamReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    po::options_description bam_desc("Headless SAM/BAM Output");
    bam_desc.add_options()("o-hl_sam", po::value<std::string>(),
        "Destination of output headless SAM/BAM file. Unset to disable the writer.");
    bam_desc.add_options()("o-hl_sam-use_m", "Whether to use CIGAR 'M' instead of '=/X' for alignment");
    bam_desc.add_options()("o-hl_sam-write_bam", "Enforce BAM instead of SAM output.");
    bam_desc.add_options()(
        "o-hl_sam-num_threads", po::value<int>()->default_value(4), "Number of threads used in BAM compression.");
    bam_desc.add_options()("o-hl_sam-compress_level", po::value<char>()->default_value('4'),
        "Compression level in BAM. Support `u` for uncompressed raw BAM output and [0-9] for underlying zlib "
        "compression.");
    desc.add(bam_desc);
}

BaseReadOutput* HeadlessBamReadOutputFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const
{
    if (vm.count("o-hl_sam") != 0U) {
        if (fasta_fetch->num_seqs() != 0) {
            BOOST_LOG_TRIVIAL(warning) << "Sequences presented in the reference file. Use SAM/BAM instead of this "
                                          "headless one for better compatibility.";
        }
        auto so = BamOptions();
        so.use_m = vm.count("o-hl_sam-use_m") > 0;
        so.write_bam = vm.count("o-hl_sam-write_bam") > 0;
        so.PG_CL = boost::algorithm::join(args, " ");
        so.hts_io_threads = vm["o-hl_sam-num_threads"].as<int>();
        so.compress_level = vm["o-hl_sam-compress_level"].as<char>();
        if (ALLOWED_COMPRESSION_LEVELS.find(so.compress_level) == std::string::npos) {
            BOOST_LOG_TRIVIAL(fatal) << "Invalid compression level: " << so.compress_level
                                     << ". Allowed values are: " << ALLOWED_COMPRESSION_LEVELS;
            abort_mpi();
        }
        return new HeadlessBamReadOutput(vm["o-hl_sam"].as<std::string>(), so);
    }
    return new DumbReadOutput();
}

HeadlessBamReadOutputFactory::~HeadlessBamReadOutputFactory() = default;
} // namespace labw::art_modern
