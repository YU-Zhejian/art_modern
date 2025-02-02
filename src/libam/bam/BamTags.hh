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
    using tag_type = std::tuple<std::string, char, int, data_type>;
    void patch(bam1_t* record) const;
    [[nodiscard]] size_t size() const;
    void add_string(const std::string& key, const std::string& value);
    void add_int_i(const std::string& key, int32_t value);

private:
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