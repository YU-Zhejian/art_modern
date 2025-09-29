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

#pragma once

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/lockfree/SimpleLFIO.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>

#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

class PwaReadOutput final: public BaseReadOutput {
public:
    DELETE_MOVE(PwaReadOutput)
    DELETE_COPY(PwaReadOutput)

    explicit PwaReadOutput(const std::string& filename, const std::vector<std::string>& args, int n_threads);
    void writeSE(const ProducerToken& token, const PairwiseAlignment& pwa) override;
    void writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    ProducerToken get_producer_token() override;

    bool require_alignment() const override;

    void close() override;

    ~PwaReadOutput() override;

private:
    SimpleLFIO lfio_;
};

class PwaReadOutputFactory final: public BaseReadOutputFactory {
public:
    DELETE_MOVE(PwaReadOutputFactory)
    DELETE_COPY(PwaReadOutputFactory)
    PwaReadOutputFactory() = default;
    ~PwaReadOutputFactory() override = default;

    [[nodiscard]] std::string name() const override { return "PWA"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    [[nodiscard]] std::shared_ptr<BaseReadOutput> create(const OutParams& params) const override;
};
} // namespace labw::art_modern
