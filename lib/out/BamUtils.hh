#pragma once
#include "PairwiseAlignment.hh"
#include "SamOptions.hh"

#include <string>

namespace labw {
namespace art_modern {

    void assert_correct_cigar(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, const uint32_t* cigar_c_arr);
    void fill_md_nm_tag(bam1_t* b, const PairwiseAlignment& pwa);

    class BamUtils {
    public:
        BamUtils(BamUtils&& other) = delete;
        BamUtils(const BamUtils&) = delete;
        BamUtils& operator=(BamUtils&&) = delete;
        BamUtils& operator=(const BamUtils&) = delete;

        explicit BamUtils(const SamOptions& sam_options);
        std::string generate_oa_tag(const PairwiseAlignment& pwa) const;

    private:
        const SamOptions& sam_options_;
    };

} // art_modern
} // labw
