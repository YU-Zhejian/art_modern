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

#include "libam_support/out/PwaReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/fs_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/algorithm/string/join.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
namespace {
    std::string preamble(const std::vector<std::string>& args)
    {
        std::ostringstream oss;
        oss << "#PWA\n";
        oss << "#ARGS: " << boost::algorithm::join(args, " ") << "\n";
        return oss.str();
    }
} // namespace

void PwaReadOutput::writeSE(const ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    auto os = std::make_unique<std::string>(pwa.serialize());
    lfio_.push(std::move(os), token);
}

void PwaReadOutput::writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    auto os1 = std::make_unique<std::string>(pwa1.serialize(1));
    lfio_.push(std::move(os1), token);
    auto os2 = std::make_unique<std::string>(pwa2.serialize(2));
    lfio_.push(std::move(os2), token);
}

PwaReadOutput::~PwaReadOutput() { PwaReadOutput::close(); }
PwaReadOutput::PwaReadOutput(const std::string& filename, const std::vector<std::string>& args, const int n_threads)
    : lfio_("PWA", filename, preamble(args))
{
    lfio_.init_queue(n_threads, 0);
    lfio_.start();
}

void PwaReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    closed_ = true;
}

bool PwaReadOutput::require_alignment() const { return true; }

ProducerToken PwaReadOutput::get_producer_token() { return lfio_.get_producer_token(); }

void PwaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description pwa_desc("PWA Output");
    pwa_desc.add_options()("o-pwa", boost::program_options::value<std::string>(),
        "Destination of output pwa file. Unset to disable the writer.");
    desc.add(pwa_desc);
}
std::shared_ptr<BaseReadOutput> PwaReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-pwa") != 0U) {
        return std::make_shared<PwaReadOutput>(
            attach_mpi_rank_to_path(params.vm["o-pwa"].as<std::string>(), mpi_rank()), params.args, params.n_threads);
    }
    throw OutputNotSpecifiedException {};
}
} // namespace labw::art_modern
