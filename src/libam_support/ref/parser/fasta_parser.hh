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
#include <string>
#include <utility>

namespace labw::art_modern {

using FastaRecord = std::pair<std::string, std::string>;

class EOFException final : std::exception {
};
class MalformedFastaException final : std::exception {
    [[nodiscard]] const char* what() const noexcept override;
};

class FastaIterator {
public:
    explicit FastaIterator(std::istream& istream);

    DELETE_COPY(FastaIterator)
    DELETE_MOVE(FastaIterator)
    ~FastaIterator() = default;

    FastaRecord next();

private:
    std::istream& _istream;
    std::size_t _lineno = 0;
    std::mutex mutex_;
    std::string staged_next_record_id_ {};
};

} // namespace labw::art_modern