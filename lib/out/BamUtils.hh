#pragma once
#include "PairwiseAlignment.hh"
#include <string>

namespace labw {
namespace art_modern {
    class SamReadOutputOptions {
    public:
        /**
         * Format version. Accepted format: `/^[0-9]+\.[0-9]+$`.
         */
        std::string HD_VN = "1.4";
        /**
         * Sorting order of alignments. Valid values: `unknown` (default), `unsorted`, `queryname` and `coordinate`.
         */
        std::string HD_SO = "unsorted";

        /**
         * Program record identifier.
         */
        std::string PG_ID = "01";
        /**
         * Program name.
         */
        std::string PG_PN = "art_modern";
        /**
         * Command line.
         */
        std::string PG_CL;
        /**
         * Program version.
         */
        std::string PG_VN = "1.0"; // "ART-" ART_VERSION "-ART_MODERN-" ART_MODERN_VERSION;

        bool use_m = false;

        /**
         * If `false`, will write SAM instead.
         */
        bool write_bam = true;
    };
    void assert_correct_cigar(
        const PairwiseAlignment& pwa, const std::vector<uint32_t>& cigar, const uint32_t* cigar_c_arr);

    class BamUtils {
    public:
        explicit BamUtils(const SamReadOutputOptions& sam_options);
        std::string generate_oa_tag(const PairwiseAlignment& pwa) const;

    private:
        const SamReadOutputOptions& sam_options_;
    };

} // art_modern
} // labw
