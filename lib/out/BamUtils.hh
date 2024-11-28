#pragma once
#include "LockFreeIO.hh"
#include "PairwiseAlignment.hh"
#include "SamOptions.hh"
#include <htslib/sam.h>

#include <string>

namespace labw::art_modern {

void assert_correct_cigar([[maybe_unused]] const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar);

class BamUtils {
public:
    BamUtils(BamUtils&& other) = delete;
    BamUtils(const BamUtils&) = delete;
    BamUtils& operator=(BamUtils&&) = delete;
    BamUtils& operator=(const BamUtils&) = delete;

    static std::string generate_oa_tag(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, int32_t nm_tag);
    static std::pair<int32_t, std::string> generate_nm_md_tag(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar);
    static bam1_t* init();
    static void write(samFile* fp, const sam_hdr_t* h, const bam1_t* b);
};

class BamTags {
public:
    using data_type = std::shared_ptr<uint8_t[]>;
    using tag_type = std::tuple<std::string, char, int, data_type>;
    void patch(bam1_t* record) const;
    size_t size() const;
    void add_string(const std::string& key, const std::string& value);
    void add_int_i(const std::string& key, int32_t value);

private:
    std::vector<tag_type> tags_;
};

class BamLFIO : public LockFreeIO<bam1_t> {
public:
    void write(bam1_t* ss) override
    {
        BamUtils::write(fp_, h_, ss);
        bam_destroy1(ss);
    }
    BamLFIO(samFile* fp, sam_hdr_t* h)
        : fp_(fp)
        , h_(h)
    {
    }

private:
    samFile* fp_;
    sam_hdr_t* h_;
};

} // namespace labw::art_modern
