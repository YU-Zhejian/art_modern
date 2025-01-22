#include "art_modern_config.h"

#include "libam/bam/BamTags.hh"

#include "libam/CExceptionsProxy.hh"

#include <htslib/sam.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {
void BamTags::patch(bam1_t* record) const
{
    for (const auto& [tag_name, tag_type, tag_len, tag_data] : tags_) {
        CExceptionsProxy::assert_numeric(bam_aux_append(
            record, tag_name.c_str(), tag_type, tag_len,
             reinterpret_cast<const uint8_t*>(tag_data->c_str())
        ),
            USED_HTSLIB_NAME, "Failed to add tag to read", false, CExceptionsProxy::EXPECTATION::ZERO);
    }
}
size_t BamTags::size() const
{
    size_t size = 0;
    for (const auto& [tag_name, tag_type, tag_len, tag_data] : tags_) {
        size += tag_len + 3;
    }
    return size;
}
void BamTags::add_string(const std::string& key, const std::string& value)
{
    const auto len = value.size();
    auto data = std::make_shared<std::string>();
    data->resize(len + 1);
    std::memcpy(data->data(), value.data(), len);
    (*data)[static_cast<std::ptrdiff_t>(len)] = 0;
    tags_.emplace_back(key, 'Z', len + 1, data);
}
void BamTags::add_int_i(const std::string& key, const int32_t value)
{
    auto data = std::make_shared<std::string>();
    data->resize(4);
    std::memcpy(data->data(), &value, 4);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), 4);
#endif
    tags_.emplace_back(key, 'i', 4, data);
}
} // namespace labw::art_modern
