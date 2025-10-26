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

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <htslib/hts.h>

#include <cstddef>
#include <istream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace labw::art_modern {

/**
 * This class has move constructor only.
 */
class InMemoryFastaFetch final : public BaseFastaFetch {
public:
    InMemoryFastaFetch(InMemoryFastaFetch&& other) noexcept
        : BaseFastaFetch(std::move(other.seq_names_), std::move(other.seq_lengths_))
        , seqs_(std::move(other.seqs_))
    {
    }
    DELETE_MOVE_ASSIGNMENT(InMemoryFastaFetch)

    DELETE_COPY(InMemoryFastaFetch)

    InMemoryFastaFetch() = default;
    ~InMemoryFastaFetch() override = default;

    InMemoryFastaFetch(const InMemoryFastaFetch& other, std::ptrdiff_t from, std::ptrdiff_t to);
    explicit InMemoryFastaFetch(const std::string& file_name);
    explicit InMemoryFastaFetch(std::istream& iss);
    explicit InMemoryFastaFetch(std::tuple<std::vector<std::string>, std::vector<std::string>> seq_map);
    InMemoryFastaFetch(std::vector<std::string>&& seq_name, std::vector<std::string>&& seq);
    std::string fetch(size_t seq_id, hts_pos_t start, hts_pos_t end) override;
    std::string fetch(size_t seq_id) override;
    void clear() override;

private:
    std::vector<std::string> seqs_;
};
} // namespace labw::art_modern
