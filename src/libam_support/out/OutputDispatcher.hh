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
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>

#include <atomic>
#include <memory>
#include <vector>

namespace labw::art_modern {

class OutputDispatcher {

public:
    DELETE_MOVE(OutputDispatcher)
    DELETE_COPY(OutputDispatcher)
    using TokenRing = std::vector<ProducerToken>;

    [[nodiscard]] bool require_alignment() const;

    OutputDispatcher() = default;
    ~OutputDispatcher();

    void add(std::shared_ptr<BaseReadOutput>&& output);
    void writeSE(const TokenRing& tokens, const PairwiseAlignment& pwa);
    void writePE(const TokenRing& tokens, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2);
    void close();
    TokenRing get_producer_tokens();

private:
    std::vector<std::shared_ptr<BaseReadOutput>> outputs_;
    std::atomic<bool> closed_ = false;
};

class OutputDispatcherFactory {
public:
    DELETE_MOVE(OutputDispatcherFactory)
    DELETE_COPY(OutputDispatcherFactory)

    OutputDispatcherFactory();

    void patch_options(boost::program_options::options_description& desc) const;
    [[nodiscard]] std::shared_ptr<OutputDispatcher> create(const OutParams& params) const;
    ~OutputDispatcherFactory();

private:
    std::vector<std::shared_ptr<BaseReadOutputFactory>> factories_;
};

} // namespace labw::art_modern