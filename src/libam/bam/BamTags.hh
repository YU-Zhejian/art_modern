#pragma once

#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {
class BamTags {
public:
    using data_type = std::shared_ptr<std::string>;
    using tag_type = std::tuple<std::string, char, int, data_type>;
    void patch(bam1_t* record) const;
    [[nodiscard]] size_t size() const;
    void add_string(const std::string& key, const std::string& value);
    void add_int_i(const std::string& key, int32_t value);

private:
    std::vector<tag_type> tags_;
};

} // namespace labw::art_modern