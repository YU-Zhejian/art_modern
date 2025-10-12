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

#include "libam_support/bam/BamLFIO.hh"
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>

#include <htslib/sam.h>

#include <memory>
#include <string>

namespace labw::art_modern {

class HeadlessBamReadOutput final : public BaseReadOutput {
public:
    DELETE_MOVE(HeadlessBamReadOutput)
    DELETE_COPY(HeadlessBamReadOutput)

    [[nodiscard]] bool require_alignment() const override;

    HeadlessBamReadOutput(const std::string& filename, const BamOptions& sam_options, std::size_t n_threads);
    void writeSE(const ProducerToken& token, const PairwiseAlignment& pwa) override;
    void writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    ProducerToken get_producer_token() override;
    void close() override;
    ~HeadlessBamReadOutput() override;

private:
    samFile* sam_file_;
    sam_hdr_t* sam_header_;
    const BamOptions sam_options_;
    BamLFIO lfio_;
};

class HeadlessBamReadOutputFactory final : public BaseReadOutputFactory {
public:
    DELETE_MOVE(HeadlessBamReadOutputFactory)
    DELETE_COPY(HeadlessBamReadOutputFactory)
    HeadlessBamReadOutputFactory() = default;

    [[nodiscard]] std::string name() const override { return "HeadlessBam"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    [[nodiscard]] std::shared_ptr<BaseReadOutput> create(const OutParams& params) const override;
    ~HeadlessBamReadOutputFactory() override;

private:
    BamOptions sam_options_;
};

} // namespace labw::art_modern
