/**
 * Copyright 2024-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "libam_support/out/HeadlessBamReadOutput.hh"

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Constants.hh"
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/bam/BamTags.hh"
#include "libam_support/bam/BamUtils.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <htslib/sam.h>

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

namespace po = boost::program_options;

namespace labw::art_modern {
HeadlessBamReadOutput::HeadlessBamReadOutput(const std::string& filename, const BamOptions& sam_options,
    const std::size_t n_threads, const std::size_t queue_size)
    : sam_file_(BamUtils::open_file(filename, sam_options))
    , sam_header_(BamUtils::init_header(sam_options))
    , sam_options_(sam_options)
    , lfio_("HeadlessBAM", sam_file_, sam_header_)
{
    CExceptionsProxy::assert_numeric(
        sam_hdr_write(sam_file_, sam_header_), USED_HTSLIB_NAME, "Failed to write SAM/BAM record");
    lfio_.init_queue(n_threads, 0, queue_size);
    lfio_.start();
}
void HeadlessBamReadOutput::writeSE(const ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    auto sam_record = BamUtils::init_uptr();
    const auto rlen = pwa.query.size();
    auto seq = pwa.query;
    auto cigar = pwa.generate_cigar_array(sam_options_.use_m);
    BamUtils::assert_correct_cigar(pwa, cigar);
    if (!pwa.is_plus_strand) {
        reverse(cigar.data(), cigar.size());
        revcomp_inplace(seq);
    }
    BamTags tags;
    if (sam_options_.with_tag_MD || sam_options_.with_tag_NM || sam_options_.with_tag_OA) {
        const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa, cigar);
        if (sam_options_.with_tag_MD) {
            tags.add_string("MD", md_tag);
        }
        if (sam_options_.with_tag_NM) {
            tags.add_int_i("NM", nm_tag);
        }
        if (sam_options_.with_tag_OA) {
            const auto& oa_tag = BamUtils::generate_oa_tag(pwa, cigar, nm_tag);
            tags.add_string("OA", oa_tag);
        }
    }

    CExceptionsProxy::assert_numeric(
        bam_set1(sam_record.get(), pwa.read_name.size(), pwa.read_name.c_str(),
            BAM_FUNMAP, // Alignment info moved to OA tag
            TID_FOR_UNMAPPED, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            0, // Alignment info moved to OA tag
            nullptr, // Alignment info moved to OA tag
            TID_FOR_UNMAPPED, // Unset for SE reads
            0, // Unset for SE reads
            0, // Unset for SE reads
            rlen, seq.c_str(), sam_options_.no_qual ? nullptr : reinterpret_cast<const char*>(pwa.qual_vec.data()),
            tags.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);
    if (!sam_options_.no_qual) {
        if (!pwa.is_plus_strand) {
            reverse(bam_get_qual(sam_record), rlen);
        }
    }
    tags.patch(sam_record.get());
    lfio_.push(std::move(sam_record), token);
}
void HeadlessBamReadOutput::writePE(
    const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    auto sam_record1 = BamUtils::init_uptr();
    auto sam_record2 = BamUtils::init_uptr();

    const auto rlen_1 = pwa1.query.size();
    const auto rlen_2 = pwa2.query.size();

    auto seq1 = pwa1.query;
    auto seq2 = pwa2.query;

    auto cigar1 = pwa1.generate_cigar_array(sam_options_.use_m);
    auto cigar2 = pwa2.generate_cigar_array(sam_options_.use_m);

    BamUtils::assert_correct_cigar(pwa1, cigar1);
    BamUtils::assert_correct_cigar(pwa2, cigar2);

    if (!pwa1.is_plus_strand) {
        reverse(cigar1.data(), cigar1.size());
        revcomp_inplace(seq1);
    }
    if (!pwa2.is_plus_strand) {
        reverse(cigar2.data(), cigar2.size());
        revcomp_inplace(seq2);
    }

    BamTags tags1;
    if (sam_options_.with_tag_MD || sam_options_.with_tag_NM || sam_options_.with_tag_OA) {
        const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa1, cigar1);
        if (sam_options_.with_tag_MD) {
            tags1.add_string("MD", md_tag);
        }
        if (sam_options_.with_tag_NM) {
            tags1.add_int_i("NM", nm_tag);
        }
        if (sam_options_.with_tag_OA) {
            const auto& oa_tag = BamUtils::generate_oa_tag(pwa1, cigar1, nm_tag);
            tags1.add_string("OA", oa_tag);
        }
    }
    BamTags tags2;
    if (sam_options_.with_tag_MD || sam_options_.with_tag_NM || sam_options_.with_tag_OA) {
        const auto& [nm_tag, md_tag] = BamUtils::generate_nm_md_tag(pwa2, cigar2);
        if (sam_options_.with_tag_MD) {
            tags2.add_string("MD", md_tag);
        }
        if (sam_options_.with_tag_NM) {
            tags2.add_int_i("NM", nm_tag);
        }
        if (sam_options_.with_tag_OA) {
            const auto& oa_tag = BamUtils::generate_oa_tag(pwa1, cigar1, nm_tag);
            tags2.add_string("OA", oa_tag);
        }
    }

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
            rlen_1, seq1.c_str(), sam_options_.no_qual ? nullptr : reinterpret_cast<const char*>(pwa1.qual_vec.data()),
            tags1.size()),
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
            rlen_2, seq2.c_str(), sam_options_.no_qual ? nullptr : reinterpret_cast<const char*>(pwa2.qual_vec.data()),
            tags2.size()),
        USED_HTSLIB_NAME, "Failed to populate SAM/BAM record", false, CExceptionsProxy::EXPECTATION::NON_NEGATIVE);

    if (!sam_options_.no_qual) {
        if (!pwa1.is_plus_strand) {
            reverse(bam_get_qual(sam_record1), rlen_1);
        } else {
            reverse(bam_get_qual(sam_record2), rlen_2);
        }
    }

    tags1.patch(sam_record1.get());
    tags2.patch(sam_record2.get());

    lfio_.push(std::move(sam_record1), token);
    lfio_.push(std::move(sam_record2), token);
}
void HeadlessBamReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    sam_hdr_destroy(sam_header_);
    closed_ = true;
}
HeadlessBamReadOutput::~HeadlessBamReadOutput() { HeadlessBamReadOutput::close(); }

bool HeadlessBamReadOutput::require_alignment() const { return true; }

ProducerToken HeadlessBamReadOutput::get_producer_token() { return lfio_.get_producer_token(); }

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
    bam_desc.add_options()("o-hl_sam-queue_size",
        boost::program_options::value<std::size_t>()->default_value(LockFreeIO<void*>::QUEUE_SIZE),
        "Size of the lock-free queue used in headless SAM/BAM output.");
    bam_desc.add_options()("o-hl_sam-without_tag_OA", "Set to disable the OA tag in headless SAM/BAM output.");
    bam_desc.add_options()("o-hl_sam-without_tag_MD", "Set to disable the MD tag in headless SAM/BAM output.");
    bam_desc.add_options()("o-hl_sam-without_tag_NM", "Set to disable the NM tag in headless SAM/BAM output.");
    bam_desc.add_options()("o-hl_sam-no_qual", "Set to disable writing quality scores in headless SAM/BAM output.");
    desc.add(bam_desc);
}

std::shared_ptr<BaseReadOutput> HeadlessBamReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-hl_sam") != 0U) {
        if (params.fasta_fetch->num_seqs() != 0) {
            BOOST_LOG_TRIVIAL(warning) << "Sequences presented in the reference file. Use SAM/BAM instead of this "
                                          "headless one for better compatibility.";
        }
        auto so = BamOptions();
        so.use_m = params.vm.count("o-hl_sam-use_m") > 0;
        so.write_bam = params.vm.count("o-hl_sam-write_bam") > 0;
        so.PG_CL = join(params.args, " ");
        so.hts_io_threads = params.vm["o-hl_sam-num_threads"].as<int>();
        so.compress_level = params.vm["o-hl_sam-compress_level"].as<char>();
        if (std::string(BamOptions::ALLOWED_COMPRESSION_LEVELS).find(so.compress_level) == std::string::npos) {
            BOOST_LOG_TRIVIAL(fatal) << "Invalid compression level: " << so.compress_level
                                     << ". Allowed values are: " << BamOptions::ALLOWED_COMPRESSION_LEVELS;
            abort_mpi();
        }
        so.with_tag_OA = params.vm.count("o-hl_sam-without_tag_OA") == 0;
        so.with_tag_MD = params.vm.count("o-hl_sam-without_tag_MD") == 0;
        so.with_tag_NM = params.vm.count("o-hl_sam-without_tag_NM") == 0;
        so.no_qual = params.vm.count("o-hl_sam-no_qual") > 0;
        so.log_("HeadlessBamReadOutput");
        return std::make_shared<HeadlessBamReadOutput>(
            attach_mpi_rank_to_path(params.vm["o-hl_sam"].as<std::string>(), mpi_rank_s()), so, params.n_threads,
            params.vm["o-hl_sam-queue_size"].as<std::size_t>());
    }
    throw OutputNotSpecifiedException {};
}
} // namespace labw::art_modern
