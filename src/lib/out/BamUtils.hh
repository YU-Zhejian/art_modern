#pragma once
#include "LockFreeIO.hh"
#include "PairwiseAlignment.hh"
#include "art_modern_dtypes.hh"
#include "SamOptions.hh"
#include <htslib/sam.h>
#include <memory>

#include <string>

namespace labw::art_modern {

void assert_correct_cigar(
    [[maybe_unused]] const PairwiseAlignment& pwa, [[maybe_unused]] const std::vector<am_cigar_t>& cigar);

struct BamDestroyer {
    void operator()(bam1_t* b) const { bam_destroy1(b); }
};

class BamUtils {
public:
    BamUtils(BamUtils&& other) = delete;
    BamUtils(const BamUtils&) = delete;
    BamUtils& operator=(BamUtils&&) = delete;
    BamUtils& operator=(const BamUtils&) = delete;
    using bam1_t_uptr = std::unique_ptr<bam1_t, BamDestroyer>;

    static std::string generate_oa_tag(
        const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar, int32_t nm_tag);
    static std::pair<int32_t, std::string> generate_nm_md_tag(
        const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar);
    static bam1_t* init();
    static bam1_t_uptr init_uptr();
    static sam_hdr_t* init_header(const SamOptions& sam_options);
    static samFile* open_file(const std::string& filename, const SamOptions& sam_options);
    static void write(samFile* fp, const sam_hdr_t* h, const bam1_t* b);
};

class BamTags {
public:
    using data_type = std::shared_ptr<uint8_t[]>;
    using tag_type = std::tuple<std::string, char, int, data_type>;
    void patch(bam1_t* record) const;
    [[nodiscard]] size_t size() const;
    void add_string(const std::string& key, const std::string& value);
    void add_int_i(const std::string& key, int32_t value);

private:
    std::vector<tag_type> tags_;
};

class BamLFIO : public LockFreeIO<BamUtils::bam1_t_uptr> {
public:
    void write(BamUtils::bam1_t_uptr ss) override { BamUtils::write(fp_, h_, ss.get()); }
    BamLFIO(samFile* fp, const sam_hdr_t* h)
        : fp_(fp)
        , h_(h)
    {
    }
    ~BamLFIO() override { stop(); };

private:
    samFile* fp_;
    const sam_hdr_t* h_;
};

} // namespace labw::art_modern
