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

#include "libam_support/utils/class_macros_utils.hh"

#include <cstddef>
#include <exception>
#include <istream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

namespace labw::art_modern {

/**
 * Indicates the end of the FASTA file.
 */
class FastaEOFException final : public std::exception { };

class MalformedFastaException final : public std::runtime_error {
public:
    explicit MalformedFastaException(const std::string& reason)
        : std::runtime_error(reason) {};
    DEFAULT_DESTRUCTOR(MalformedFastaException)
    DEFAULT_COPY(MalformedFastaException)
    DEFAULT_MOVE(MalformedFastaException)
};

class FastaIterator {
public:
    using FastaRecord = std::pair<std::string, std::string>;
    explicit FastaIterator(std::istream& istream);

    DELETE_COPY(FastaIterator)
    DELETE_MOVE(FastaIterator)
    DEFAULT_DESTRUCTOR(FastaIterator)

    FastaRecord next();

private:
    std::istream& _istream;
    std::size_t _lineno = 0;
    std::mutex mutex_;
    std::string staged_next_record_id_;
};

} // namespace labw::art_modern
