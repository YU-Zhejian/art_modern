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
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>

#include <atomic>
#include <exception>
#include <memory>
#include <string>

namespace labw::art_modern {

class OutputNotSpecifiedException final: public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override { return "Output file not specified"; }
};

class BaseReadOutput {
public:
    DELETE_MOVE(BaseReadOutput)
    DELETE_COPY(BaseReadOutput)

    BaseReadOutput() = default;
    virtual void writeSE(const ProducerToken& token, const PairwiseAlignment& pwa) = 0;
    virtual void writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) = 0;
    virtual void close() = 0;
    [[nodiscard]] virtual bool require_alignment() const = 0;
    virtual ~BaseReadOutput() = default;
    virtual ProducerToken get_producer_token() = 0;

protected:
    std::atomic<bool> closed_ = false;
};

class BaseReadOutputFactory {
public:
    DELETE_MOVE(BaseReadOutputFactory)
    DELETE_COPY(BaseReadOutputFactory)
    BaseReadOutputFactory() = default;

    [[nodiscard]] virtual std::string name() const = 0;
    virtual void patch_options(boost::program_options::options_description& desc) const = 0;
    /**
     * @throw OutputNotSpecifiedException if output file is not specified.
     */
    [[nodiscard]] virtual std::shared_ptr<BaseReadOutput> create(const OutParams& params) const = 0;
    virtual ~BaseReadOutputFactory() = default;
};

} // namespace labw::art_modern
