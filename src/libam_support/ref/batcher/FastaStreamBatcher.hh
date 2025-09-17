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
#include "libam_support/ref/parser/fasta_parser.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <istream>
#include <mutex>

namespace labw::art_modern {
class FastaStreamBatcher final {
public:
    /**
     *
     * @return A {@link InMemoryFastaFetch} of less or equal than {@link batch_size} items.
     * Empty {@link InMemoryFastaFetch} if the stream is exhausted.
     */
    InMemoryFastaFetch fetch();

    FastaStreamBatcher(std::size_t batch_size, std::istream& stream);

    DELETE_COPY(FastaStreamBatcher)
    DELETE_MOVE(FastaStreamBatcher)
    ~FastaStreamBatcher() = default;

private:
    std::size_t batch_size_;
    FastaIterator fasta_iterator_;
    std::mutex mutex_;
};

} // namespace labw::art_modern