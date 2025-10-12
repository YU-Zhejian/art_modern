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

#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <memory>
#include <mutex>

namespace labw::art_modern {

class InMemoryFastaBatcher {
public:
    InMemoryFastaBatcher(std::size_t batch_size, const std::shared_ptr<InMemoryFastaFetch>& stream);
    ~InMemoryFastaBatcher() = default;

    DELETE_COPY(InMemoryFastaBatcher)
    DELETE_MOVE(InMemoryFastaBatcher)
    /**
     *
     * @return A {@link InMemoryFastaFetch @endlink} of less or equal than {@link batch_size @endlink} items.
     * Empty {@link InMemoryFastaFetch @endlink} if the stream is exhausted.
     */
    std::shared_ptr<InMemoryFastaFetch> fetch();

private:
    std::size_t batch_size_;
    std::size_t current_index_ { 0 };
    const std::shared_ptr<InMemoryFastaFetch>& stream_;
    std::mutex mutex_;
};

} // namespace labw::art_modern
