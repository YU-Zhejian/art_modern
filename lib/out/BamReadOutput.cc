#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>
#include <utility>

#include "BamReadOutput.hh"
#include "CExceptionsProxy.hh"
#include "DumbReadOutput.hh"
#include "art_modern_config.h" // For USED_HTSLIB_NAME
#include "utils/mpi_utils.hh"
#include "utils/seq_utils.hh"

namespace po = boost::program_options;

namespace labw::art_modern {

void BamReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }

    const int tid = CExceptionsProxy::assert_numeric(sam_hdr_name2tid(sam_header_, pwa.contig_name.c_str()),
        USED_HTSLIB_NAME, "Failed to fetch TID for contig '" + pwa.contig_name + "'", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    auto sam_record = BamUtils::init();
    const auto rlen = static_cast<long>(pwa.query.size());
    auto cigar = pwa.generate_cigar_array(sam_options_.use_m);
    assert_correct_cigar(pwa, cigar);

    const auto& seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);
    auto qual = pwa.qual;
    const hts_pos_t pos = pwa.pos_on_contig;
    if (!pwa.is_plus_strand) {
        std::reverse(qual.begin(), qual.end());
        std::reverse(cigar.begin(), cigar.end());
    }

    const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa, cigar);
    BamTags tags;
    tags.add_string("MD", md_tag);
    tags.add_int_i("NM", nm_tag);

    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record, pwa.read_name.length(), pwa.read_name.c_str(), pwa.is_plus_strand ? 0 : BAM_FREVERSE, tid,
            pos, MAPQ_MAX, cigar.size(), cigar.data(),
            0, // Unset for SE reads
            0, // Unset for SE reads
            0, // Unset for SE reads
            rlen, seq.c_str(), qual.c_str(), tags.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    tags.patch(sam_record);

    std::unique_lock rhs_lk(mutex_);
    BamUtils::write(sam_file_, sam_header_, sam_record);
    bam_destroy1(sam_record);
}

void BamReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }

    const int tid = CExceptionsProxy::assert_numeric(sam_hdr_name2tid(sam_header_, pwa1.contig_name.c_str()),
        USED_HTSLIB_NAME, "Failed to fetch TID for contig '" + pwa1.contig_name + "'", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    const auto rlen = static_cast<long>(pwa1.query.size());

    auto cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
    auto cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

    assert_correct_cigar(pwa1, cigar1);
    assert_correct_cigar(pwa2, cigar2);

    const auto& [nm_tag1, md_tag1] = BamUtils::generate_nm_md_tag(pwa1, cigar1);
    const auto& [nm_tag2, md_tag2] = BamUtils::generate_nm_md_tag(pwa2, cigar2);

    BamTags tags1;
    tags1.add_string("MD", md_tag1);
    tags1.add_int_i("NM", nm_tag1);

    BamTags tags2;
    tags2.add_string("MD", md_tag2);
    tags2.add_int_i("NM", nm_tag2);

    const auto& seq1 = pwa1.is_plus_strand ? pwa1.query : revcomp(pwa1.query);
    const auto& seq2 = pwa2.is_plus_strand ? pwa2.query : revcomp(pwa2.query);

    const hts_pos_t pos1 = pwa1.pos_on_contig;
    const hts_pos_t pos2 = pwa2.pos_on_contig;

    auto qual1 = pwa1.qual;
    auto qual2 = pwa2.qual;

    uint16_t flag1 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
    uint16_t flag2 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;

    if (pwa1.is_plus_strand) {
        flag1 |= BAM_FMREVERSE;
        flag2 |= BAM_FREVERSE;
        std::reverse(qual2.begin(), qual2.end());
        std::reverse(cigar2.begin(), cigar2.end());
    } else {
        flag1 |= BAM_FREVERSE;
        flag2 |= BAM_FMREVERSE;
        std::reverse(qual1.begin(), qual1.end());
        std::reverse(cigar1.begin(), cigar1.end());
    }

    const hts_pos_t isize1 = pos2 > pos1 ? pos2 + rlen - pos1 : -(pos1 + rlen - pos2);
    const hts_pos_t isize2 = -isize1;

    auto sam_record1 = BamUtils::init();
    auto sam_record2 = BamUtils::init();

    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record1, pwa1.read_name.length(), pwa1.read_name.c_str(), flag1, tid, pos1, MAPQ_MAX,
            cigar1.size(), cigar1.data(), tid, pos2, isize1, rlen, seq1.c_str(), qual1.c_str(), tags1.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record2, pwa2.read_name.length(), pwa2.read_name.c_str(), flag2, tid, pos2, MAPQ_MAX,
            cigar2.size(), cigar2.data(), tid, pos1, isize2, rlen, seq2.c_str(), qual2.c_str(), tags2.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    tags1.patch(sam_record1);
    tags2.patch(sam_record2);

    std::scoped_lock rhs_lk(mutex_);
    BamUtils::write(sam_file_, sam_header_, sam_record1);
    BamUtils::write(sam_file_, sam_header_, sam_record2);
    bam_destroy1(sam_record1);
    bam_destroy1(sam_record2);
}
BamReadOutput::~BamReadOutput() { BamReadOutput::close(); }
BamReadOutput::BamReadOutput(const std::string& filename, const BaseFastaFetch* fasta_fetch, SamOptions sam_options)
    : sam_options_(std::move(sam_options))
    , filename(filename)
{
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

    fasta_fetch->update_sam_header(sam_header_);
    CExceptionsProxy::assert_numeric(
        sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
}

void BamReadOutput::close()
{

    if (is_closed_) {
        return;
    }
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' closed.";
    sam_close(sam_file_);
    sam_hdr_destroy(sam_header_);
    is_closed_ = true;
}
void BamReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    po::options_description bam_desc("SAM/BAM Output");
    bam_desc.add_options()(
        "o-sam", po::value<std::string>(), "Destination of output SAM/BAM file. Unset to disable the writer.");
    bam_desc.add_options()("o-sam-use_m", "Whether to use CIGAR 'M' instead of '=/X' for alignment");
    bam_desc.add_options()("o-sam-write_bam", "Enforce BAM instead of SAM output.");
    desc.add(bam_desc);
}
BaseReadOutput* BamReadOutputFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const
{
    if (vm.count("o-sam")) {
        if (fasta_fetch->num_seqs() == 0) {
            BOOST_LOG_TRIVIAL(error) << "No sequences in the reference file. If you used " << INPUT_FILE_PARSER_STREAM
                                     << " input parser, you should use headless SAM/BAM instead of this one.";
            abort_mpi();
        }
        auto so = SamOptions();
        so.use_m = vm.count("o-sam-use_m") > 0;
        so.write_bam = vm.count("o-sam-write_bam") > 0;
        so.PG_CL = boost::algorithm::join(args, " ");
        return new BamReadOutput(vm["o-sam"].as<std::string>(), fasta_fetch, so);
    }
    return new DumbReadOutput();
}

BamReadOutputFactory::~BamReadOutputFactory() = default;
}
