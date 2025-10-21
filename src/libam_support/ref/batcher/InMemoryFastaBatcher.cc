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

#include "libam_support/ref/batcher/InMemoryFastaBatcher.hh"

#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <mutex>

namespace labw::art_modern {
std::shared_ptr<InMemoryFastaFetch> InMemoryFastaBatcher::fetch()
{
    const std::scoped_lock lock(mutex_);
    const auto from = current_index_;
    const auto to = std::min(current_index_ + batch_size_, stream_->num_seqs());
    auto fetch
        = std::make_shared<InMemoryFastaFetch>(*stream_, static_cast<ptrdiff_t>(from), static_cast<ptrdiff_t>(to));
    if (!fetch->empty()) {
        BOOST_LOG_TRIVIAL(info) << "FASTA Read batch " << fetch->seq_name(0) << " to "
                                << fetch->seq_name(fetch->num_seqs() - 1) << " (" << fetch->num_seqs() << "） created";
    }
    current_index_ = to;
    return fetch;
}

InMemoryFastaBatcher::InMemoryFastaBatcher(
    const std::size_t batch_size, const std::shared_ptr<InMemoryFastaFetch>& stream)
    : batch_size_(batch_size)
    , stream_(stream)
{
}
} // namespace labw::art_modern
