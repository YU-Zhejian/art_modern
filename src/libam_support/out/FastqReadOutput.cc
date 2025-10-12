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

#include "libam_support/out/FastqReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <fmt/format.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <cstdio>
#include <memory>
#include <string>

namespace labw::art_modern {

namespace {
    std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa)
    {
        return std::make_unique<std::string>(fmt::format("@{}\n{}\n+\n{}\n", pwa.read_name, pwa.query, pwa.qual_str));
    }

    std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        return std::make_unique<std::string>(fmt::format("@{}/1\n{}\n+\n{}\n@{}/2\n{}\n+\n{}\n", pwa1.read_name,
            pwa1.query, pwa1.qual_str, pwa2.read_name, pwa2.query, pwa2.qual_str));
    }

} // namespace

void FastqReadOutput::writeSE(const ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fastq(pwa), token);
}

void FastqReadOutput::writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fastq(pwa1, pwa2), token);
}

FastqReadOutput::~FastqReadOutput() { FastqReadOutput::close(); }
FastqReadOutput::FastqReadOutput(const std::string& filename, const std::size_t n_threads)
    : lfio_("FASTQ", filename)
{
    lfio_.init_queue(n_threads, 0);
    lfio_.start();
}

void FastqReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    closed_ = true;
}

bool FastqReadOutput::require_alignment() const { return false; }

ProducerToken FastqReadOutput::get_producer_token() { return lfio_.get_producer_token(); }

void FastqReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description fastq_desc("FASTQ Output");
    fastq_desc.add_options()("o-fastq", boost::program_options::value<std::string>(),
        "Destination of output FASTQ file. Unset to disable the writer.");
    desc.add(fastq_desc);
}
std::shared_ptr<BaseReadOutput> FastqReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-fastq") != 0) {
        return std::make_shared<FastqReadOutput>(
            attach_mpi_rank_to_path(params.vm["o-fastq"].as<std::string>(), mpi_rank()), params.n_threads);
    }
    throw OutputNotSpecifiedException {};
}

} // namespace labw::art_modern
