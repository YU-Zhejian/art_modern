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

#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ds/SkipLoaderSettings.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <istream>
#include <memory>
#include <mutex>
#include <utility>

namespace labw::art_modern {

class Pbsim3TranscriptBatcher final {
public:
    DELETE_COPY(Pbsim3TranscriptBatcher)
    DELETE_MOVE(Pbsim3TranscriptBatcher)
    DEFAULT_DESTRUCTOR(Pbsim3TranscriptBatcher)

    explicit Pbsim3TranscriptBatcher(std::size_t batch_size, std::istream& istream, const SkipLoaderSettings& sls);
    std::pair<std::shared_ptr<InMemoryFastaFetch>, std::shared_ptr<CoverageInfo>> fetch();

private:
    std::size_t batch_size_;
    std::istream& istream_;
    std::mutex mutex_;
    const SkipLoaderSettings& sls_;
};

} // namespace labw::art_modern
