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

#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace labw::art_modern {
class BamTags {
public:
    explicit BamTags(int est_num_tags = 16);
    /**
     * The type of the data stored in the tag.
     *
     * FIXME: This still seems stupid.
     */
    using data_type = std::shared_ptr<std::string>;
    using key_type = std::string;
    using tag_type = std::tuple<key_type, char, int, data_type>;
    void patch(bam1_t* record) const;
    [[nodiscard]] size_t size() const;
    void add_string(const key_type& key, const std::string& value);
    void add_int_i(const key_type& key, std::int32_t value);

private:
    constexpr static std::size_t TAG_SIZE = 4;
    std::vector<tag_type> tags_;
    /**
     * Accumulated size.
     */
    std::size_t size_ = 0;
    /**
     * Length of tag name (2) and tag type (1).
     */
    static constexpr std::size_t size_of_tag_name_and_type = 3;
};

} // namespace labw::art_modern