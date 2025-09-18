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

#include "libam_support/out/FastaReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"

#include <fmt/format.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <memory>
#include <string>

namespace labw::art_modern {

namespace {
    std::unique_ptr<std::string> format_fasta(const PairwiseAlignment& pwa)
    {
        return std::make_unique<std::string>(fmt::format(">{}\n{}\n", pwa.read_name, pwa.query));
    }

    std::unique_ptr<std::string> format_fasta(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        return std::make_unique<std::string>(
            fmt::format(">{}\n{}\n>{}\n{}\n", pwa1.read_name, pwa1.query, pwa2.read_name, pwa2.query));
    }

} // namespace

void FastaReadOutput::writeSE(const ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fasta(pwa), token);
}

void FastaReadOutput::writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fasta(pwa1, pwa2), token);
}

FastaReadOutput::~FastaReadOutput() { FastaReadOutput::close(); }
FastaReadOutput::FastaReadOutput(const std::string& filename, const int n_threads)
    : lfio_("FASTA", filename)
{
    lfio_.init_queue(n_threads, 0);
    lfio_.start();
}

void FastaReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    closed_ = true;
}

bool FastaReadOutput::require_alignment() const { return false; }

ProducerToken FastaReadOutput::get_producer_token() { return lfio_.get_producer_token(); }

void FastaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description fasta_desc("FASTA Output");
    fasta_desc.add_options()("o-fasta", boost::program_options::value<std::string>(),
        "Destination of output FASTA file. Unset to disable the writer.");
    desc.add(fasta_desc);
}
std::shared_ptr<BaseReadOutput> FastaReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-fasta") != 0) {
        return std::make_shared<FastaReadOutput>(params.vm["o-fasta"].as<std::string>(), params.n_threads);
    }
    throw OutputNotSpecifiedException {};
}
} // namespace labw::art_modern
