#include "art_modern_config.h"

#include "libam/bam/BamTags.hh"

#include "libam/CExceptionsProxy.hh"

#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>

namespace labw::art_modern {
void BamTags::patch(bam1_t* record) const
{
    for (const auto& [tag_name, tag_type, tag_len, tag_data] : tags_) {
        CExceptionsProxy::assert_numeric(bam_aux_append(record, tag_name.c_str(), tag_type, tag_len,
                                             reinterpret_cast<const uint8_t*>(tag_data->c_str())),
            USED_HTSLIB_NAME, "Failed to add tag to read", false, CExceptionsProxy::EXPECTATION::ZERO);
    }
}
size_t BamTags::size() const { return size_; }
void BamTags::add_string(const std::string& key, const std::string& value)
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
void BamTags::add_int_i(const std::string& key, const int32_t value)
{
    const size_t tag_size = 4;
    auto data = std::make_shared<std::string>();
    data->resize(tag_size);
    std::memcpy(data->data(), &value, tag_size);
#ifdef HTS_LITTLE_ENDIAN
#else
    reverse(data.get(), tag_size);
#endif
    tags_.emplace_back(key, 'i', tag_size, data);
    size_ += tag_size + size_of_tag_name_and_type;
}
} // namespace labw::art_modern
