#pragma once
#include "BamReadOutput.hh"
#include "PairwiseAlignment.hh"
#include <string>

namespace labw {
namespace art_modern {

    class BamUtils {
    public:
        explicit BamUtils(const SamReadOutputOptions& sam_options);
        std::string generate_oa_tag(const PairwiseAlignment& pwa) const;

    private:
        const SamReadOutputOptions& sam_options_;
    };

} // art_modern
} // labw
