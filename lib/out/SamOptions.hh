#pragma once
#include "art_modern_config.h"
#include "art_modern_constants.hh"
#include <string>

namespace labw {
namespace art_modern {
    struct SamOptions {
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
        std::string PG_VN = std::string("ART-") + ART_VERSION + "-ART_MODERN-" + ART_MODERN_VERSION;

        /**
         * Use `M` instead of `X` or `=` for matching.
         */
        bool use_m = false;

        /**
         * If `false`, will write SAM instead.
         */
        bool write_bam = true;
    };
}
}