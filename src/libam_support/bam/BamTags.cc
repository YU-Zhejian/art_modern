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

#include "art_modern_config.h"

#include "libam_support/bam/BamTags.hh"

#include "libam_support/CExceptionsProxy.hh"

#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>

namespace labw::art_modern {
BamTags::BamTags(const int est_num_tags) { tags_.reserve(est_num_tags); }
void BamTags::patch(bam1_t* record) const
{
    for (const auto& [tag_name, tag_type, tag_len, tag_data] : tags_) {
        std::uint8_t* data = std::calloc(tag_data->size(), sizeof(std::uint8_t));
        std::memcpy(data, tag_data->data(), tag_data->size());
        CExceptionsProxy::assert_numeric(bam_aux_append(record, tag_name.c_str(), tag_type, tag_len, data),
            USED_HTSLIB_NAME, "Failed to add tag to read", false, CExceptionsProxy::EXPECTATION::ZERO);
        std::free(data);
    }
}
size_t BamTags::size() const { return size_; }
void BamTags::add_string(const key_type& key, const std::string& value)
{
    const auto len = value.size();
    const size_t tag_size = len + 1;
    auto data = std::make_shared<std::string>();
    data->resize(tag_size);
    std::memcpy(data->data(), value.data(), len);
    (*data)[static_cast<std::ptrdiff_t>(len)] = 0;
    tags_.emplace_back(key, 'Z', tag_size, data);
    size_ += tag_size + size_of_tag_name_and_type;
}
void BamTags::add_int_i(const key_type& key, const std::int32_t value)
{
    auto data = std::make_shared<std::string>();
    data->resize(TAG_SIZE);
    std::memcpy(data->data(), &value, TAG_SIZE);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), TAG_SIZE);
#endif
    tags_.emplace_back(key, 'i', TAG_SIZE, data);
    size_ += TAG_SIZE + size_of_tag_name_and_type;
}
} // namespace labw::art_modern
