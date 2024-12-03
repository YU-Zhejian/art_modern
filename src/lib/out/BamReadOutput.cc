#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>

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
    const auto& cigar = pwa.generate_cigar_array(sam_options_.use_m);
    assert_correct_cigar(pwa, cigar);

    const auto& seq = pwa.is_plus_strand ? pwa.query : revcomp(pwa.query);
    const hts_pos_t pos = pwa.pos_on_contig;

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
            rlen, seq.c_str(), pwa.qual.c_str(), tags.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    if (!pwa.is_plus_strand) {
        reverse(bam_get_qual(sam_record), rlen);
        reverse(bam_get_cigar(sam_record), sam_record->core.n_cigar);
    }
    tags.patch(sam_record);

    lfio_.push(sam_record);
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

    const auto& cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
    const auto& cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

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

    uint16_t flag1 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
    uint16_t flag2 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;

    if (pwa1.is_plus_strand) {
        flag1 |= BAM_FMREVERSE;
        flag2 |= BAM_FREVERSE;
    } else {
        flag1 |= BAM_FREVERSE;
        flag2 |= BAM_FMREVERSE;
    }

    const hts_pos_t isize1 = pos2 > pos1 ? pos2 + rlen - pos1 : -(pos1 + rlen - pos2);
    const hts_pos_t isize2 = -isize1;

    auto sam_record1 = BamUtils::init();
    auto sam_record2 = BamUtils::init();

    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record1, pwa1.read_name.length(), pwa1.read_name.c_str(), flag1, tid, pos1, MAPQ_MAX,
            cigar1.size(), cigar1.data(), tid, pos2, isize1, rlen, seq1.c_str(), pwa1.qual.c_str(), tags1.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record2, pwa2.read_name.length(), pwa2.read_name.c_str(), flag2, tid, pos2, MAPQ_MAX,
            cigar2.size(), cigar2.data(), tid, pos1, isize2, rlen, seq2.c_str(), pwa2.qual.c_str(), tags2.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    if (pwa1.is_plus_strand) {
        reverse(bam_get_qual(sam_record2), rlen);
        reverse(bam_get_cigar(sam_record2), sam_record2->core.n_cigar);
    } else {
        reverse(bam_get_qual(sam_record1), rlen);
        reverse(bam_get_cigar(sam_record1), sam_record1->core.n_cigar);
    }

    tags1.patch(sam_record1);
    tags2.patch(sam_record2);
    lfio_.push(sam_record1);
    lfio_.push(sam_record2);
}
BamReadOutput::~BamReadOutput() { BamReadOutput::close(); }
BamReadOutput::BamReadOutput(
    const std::string& filename, const BaseFastaFetch* fasta_fetch, const SamOptions& sam_options)
    : BaseFileReadOutput(filename)
    , sam_file_(BamUtils::open_file(filename, sam_options))
    , sam_header_(BamUtils::init_header(sam_options))
    , sam_options_(sam_options)
    , lfio_(sam_file_, sam_header_)
{
    fasta_fetch->update_sam_header(sam_header_);
    CExceptionsProxy::assert_numeric(
        sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
    lfio_.start();
}

void BamReadOutput::close()
{
    if (is_closed_) {
        return;
    }
    lfio_.stop();
    sam_close(sam_file_);
    sam_hdr_destroy(sam_header_);
    BaseFileReadOutput::close();
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
const std::string BamReadOutputFactory::name() const { return "BAM"; }

BamReadOutputFactory::~BamReadOutputFactory() = default;
}
