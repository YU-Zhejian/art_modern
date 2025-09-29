/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "art_modern_config.h" // For USED_HTSLIB_NAME

#include "libam_support/out/BamReadOutput.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Constants.hh"
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/bam/BamTags.hh"
#include "libam_support/bam/BamUtils.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp> // NOLINT
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdint>
#include <memory>
#include <string>
#include <utility>

namespace po = boost::program_options;

namespace labw::art_modern {

void BamReadOutput::writeSE(const ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    const int tid = CExceptionsProxy::assert_numeric(sam_hdr_name2tid(sam_header_, pwa.contig_name.c_str()),
        USED_HTSLIB_NAME, "Failed to fetch TID for contig '" + pwa.contig_name + "'", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    auto sam_record = BamUtils::init_uptr();
    const auto rlen = pwa.query.size();
    const auto& cigar = pwa.generate_cigar_array(sam_options_.use_m);
    BamUtils::assert_correct_cigar(pwa, cigar);

    auto seq = pwa.query;
    if (!pwa.is_plus_strand) {
        revcomp_inplace(seq);
    }
    const hts_pos_t pos = pwa.pos_on_contig;

    const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa, cigar);
    BamTags tags;
    tags.add_string("MD", md_tag);
    tags.add_int_i("NM", nm_tag);

    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record.get(), pwa.read_name.length(), pwa.read_name.c_str(), pwa.is_plus_strand ? 0 : BAM_FREVERSE,
                 tid, pos, MAPQ_MAX, cigar.size(), cigar.data(),
                 0, // Unset for SE reads
            0, // Unset for SE reads
            0, // Unset for SE reads
            rlen, seq.c_str(), reinterpret_cast<const char *>(pwa.qual_vec.data()), tags.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    if (!pwa.is_plus_strand) {
        reverse(bam_get_qual(sam_record.get()), rlen);
        reverse(bam_get_cigar(sam_record.get()), sam_record->core.n_cigar);
    }
    tags.patch(sam_record.get());

    lfio_.push(std::move(sam_record), token);
}

void BamReadOutput::writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    const int tid = CExceptionsProxy::assert_numeric(sam_hdr_name2tid(sam_header_, pwa1.contig_name.c_str()),
        USED_HTSLIB_NAME, "Failed to fetch TID for contig '" + pwa1.contig_name + "'", false,
        CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    const auto rlen = pwa1.query.size();

    const auto& cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
    const auto& cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

    BamUtils::assert_correct_cigar(pwa1, cigar1);
    BamUtils::assert_correct_cigar(pwa2, cigar2);

    const auto& [nm_tag1, md_tag1] = BamUtils::generate_nm_md_tag(pwa1, cigar1);
    const auto& [nm_tag2, md_tag2] = BamUtils::generate_nm_md_tag(pwa2, cigar2);

    BamTags tags1;
    tags1.add_string("MD", md_tag1);
    tags1.add_int_i("NM", nm_tag1);

    BamTags tags2;
    tags2.add_string("MD", md_tag2);
    tags2.add_int_i("NM", nm_tag2);

    auto seq1 = pwa1.query;
    auto seq2 = pwa2.query;

    const hts_pos_t pos1 = pwa1.pos_on_contig;
    const hts_pos_t pos2 = pwa2.pos_on_contig;

    uint16_t flag1 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD1;
    uint16_t flag2 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2;

    if (pwa1.is_plus_strand) {
        flag1 |= BAM_FMREVERSE;
        flag2 |= BAM_FREVERSE;
        revcomp_inplace(seq2);
    } else {
        flag1 |= BAM_FREVERSE;
        flag2 |= BAM_FMREVERSE;
        revcomp_inplace(seq1);
    }

    const hts_pos_t isize1
        = pos2 > pos1 ? static_cast<hts_pos_t>(pos2 + rlen - pos1) : -static_cast<hts_pos_t>(pos1 + rlen - pos2);
    const hts_pos_t isize2 = -isize1;

    auto sam_record1 = BamUtils::init_uptr();
    auto sam_record2 = BamUtils::init_uptr();

    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record1.get(), pwa1.read_name.length(), pwa1.read_name.c_str(), flag1, tid, pos1, MAPQ_MAX,
                 cigar1.size(), cigar1.data(), tid, pos2, isize1, rlen, seq1.c_str(),
                 reinterpret_cast<const char *>(pwa1.qual_vec.data()), tags1.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record2.get(), pwa2.read_name.length(), pwa2.read_name.c_str(), flag2, tid, pos2, MAPQ_MAX,
                 cigar2.size(), cigar2.data(), tid, pos1, isize2, rlen, seq2.c_str(), reinterpret_cast<const char *>(pwa2.qual_vec.data()), tags2.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    if (!pwa1.is_plus_strand) {
        reverse(bam_get_qual(sam_record1.get()), rlen);
        reverse(bam_get_cigar(sam_record1.get()), sam_record1->core.n_cigar);
    } else {
        reverse(bam_get_qual(sam_record2.get()), rlen);
        reverse(bam_get_cigar(sam_record2.get()), sam_record2->core.n_cigar);
    }

    tags1.patch(sam_record1.get());
    tags2.patch(sam_record2.get());

    lfio_.push(std::move(sam_record1), token);
    lfio_.push(std::move(sam_record2), token);
}
BamReadOutput::~BamReadOutput() { BamReadOutput::close(); }
BamReadOutput::BamReadOutput(const std::string& filename, const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
    const BamOptions& sam_options, const int n_threads)
    : sam_file_(BamUtils::open_file(filename, sam_options))
    , sam_header_(BamUtils::init_header(sam_options))
    , sam_options_(sam_options)
    , lfio_("BAM", sam_file_, sam_header_)
{
    fasta_fetch->update_sam_header(sam_header_);
    CExceptionsProxy::assert_numeric(
        sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
    lfio_.init_queue(n_threads, 0);
    lfio_.start();
}

void BamReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    sam_hdr_destroy(sam_header_);
    closed_ = true;
}

bool BamReadOutput::require_alignment() const { return true; }

ProducerToken BamReadOutput::get_producer_token() { return lfio_.get_producer_token(); }

void BamReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    po::options_description bam_desc("SAM/BAM Output");
    bam_desc.add_options()(
        "o-sam", po::value<std::string>(), "Destination of output SAM/BAM file. Unset to disable the writer.");
    bam_desc.add_options()("o-sam-use_m", "Whether to use CIGAR 'M' instead of '=/X' for alignment");
    bam_desc.add_options()("o-sam-write_bam", "Enforce BAM instead of SAM output.");
    bam_desc.add_options()(
        "o-sam-num_threads", po::value<int>()->default_value(4), "Number of threads used in BAM compression.");
    bam_desc.add_options()("o-sam-compress_level", po::value<char>()->default_value('4'),
        "Compression level in BAM. Support `u` for uncompressed raw BAM output and [0-9] for underlying zlib "
        "compression.");
    desc.add(bam_desc);
}
std::shared_ptr<BaseReadOutput> BamReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-sam") != 0U) {
        if (params.fasta_fetch->num_seqs() == 0) {
            BOOST_LOG_TRIVIAL(error) << "No sequences in the reference file. If you used " << INPUT_FILE_PARSER_STREAM
                                     << " input parser, you should use headless SAM/BAM instead of this one.";
            abort_mpi();
        }
        auto so = BamOptions();
        so.use_m = params.vm.count("o-sam-use_m") > 0;
        so.write_bam = params.vm.count("o-sam-write_bam") > 0;
        so.PG_CL = boost::algorithm::join(params.args, " ");
        so.hts_io_threads = params.vm["o-sam-num_threads"].as<int>();
        so.compress_level = params.vm["o-sam-compress_level"].as<char>();
        if (ALLOWED_COMPRESSION_LEVELS.find(so.compress_level) == std::string::npos) {
            BOOST_LOG_TRIVIAL(fatal) << "Invalid compression level: " << so.compress_level
                                     << ". Allowed values are: " << ALLOWED_COMPRESSION_LEVELS;
            abort_mpi();
        }
        return std::make_shared<BamReadOutput>(
            params.vm["o-sam"].as<std::string>(), params.fasta_fetch, so, params.n_threads);
    }
    throw OutputNotSpecifiedException {};
}
const std::string BamReadOutputFactory::name() const { return "BAM"; }

BamReadOutputFactory::~BamReadOutputFactory() = default;
} // namespace labw::art_modern
