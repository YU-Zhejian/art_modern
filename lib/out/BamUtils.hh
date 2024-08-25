#pragma once
#include "PairwiseAlignment.hh"
#include "SamOptions.hh"

#include <string>

namespace labw {
namespace art_modern {

    void assert_correct_cigar(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, const uint32_t* cigar_c_arr);

    class BamUtils {
    public:
        explicit BamUtils(const SamOptions& sam_options);
        std::string generate_oa_tag(const PairwiseAlignment& pwa) const;

    private:
        const SamOptions& sam_options_;
    };

} // art_modern
} // labw
