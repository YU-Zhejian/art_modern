#pragma once
#include "PairwiseAlignment.hh"
#include "SamOptions.hh"
#include <htslib/sam.h>

#include <memory>
#include <string>

namespace labw::art_modern {

void assert_correct_cigar(const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar);

class BamUtils {
public:
    BamUtils(BamUtils&& other) = delete;
    BamUtils(const BamUtils&) = delete;
    BamUtils& operator=(BamUtils&&) = delete;
    BamUtils& operator=(const BamUtils&) = delete;

    explicit BamUtils(const SamOptions& sam_options);
    static std::string generate_oa_tag(const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar);
    static std::string generate_oa_tag(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, int32_t nm_tag);
    static std::pair<int32_t, std::string> generate_nm_md_tag(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar);
    static bam1_t* init();
    static void write(samFile* fp, const sam_hdr_t* h, const bam1_t* b);

private:
    const SamOptions& sam_options_;
};

class BamTags {
public:
    using data_type = std::shared_ptr<uint8_t[]>;
    using tag_type = std::tuple<std::string, char, int, data_type>;
    void patch(bam1_t* record) const;
    size_t size() const;
    void add_string(const std::string& key, const std::string& value);
    void add_int_c(const std::string& key, int8_t value);
    void add_int_C(const std::string& key, uint8_t value);
    void add_int_s(const std::string& key, int16_t value);
    void add_int_S(const std::string& key, uint16_t value);
    void add_int_i(const std::string& key, int32_t value);
    void add_int_I(const std::string& key, uint32_t value);
    // void add_float(const std::string& key, float value);
private:
    std::vector<tag_type> tags_;
};

} // namespace labw::art_modern
