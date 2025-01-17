#pragma once

#include "libam/Dtypes.hh"
#include "libam/bam/BamOptions.hh"
#include "libam/bam/BamTypes.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/utils/class_macros_utils.hh"

#include <htslib/sam.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

class BamUtils {
public:
    DELETE_MOVE(BamUtils)
    DELETE_COPY(BamUtils)
    ~BamUtils() = default;
    BamUtils() = delete;

    static void assert_correct_cigar(
        [[maybe_unused]] const PairwiseAlignment& pwa, [[maybe_unused]] const std::vector<am_cigar_t>& cigar);
    static std::string generate_oa_tag(
        const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar, int32_t nm_tag);
    static std::pair<int32_t, std::string> generate_nm_md_tag(
        const PairwiseAlignment& pwa, const std::vector<am_cigar_t>& cigar);
    static bam1_t_uptr init_uptr();
    static sam_hdr_t* init_header(const BamOptions& sam_options);
    static samFile* open_file(const std::string& filename, const BamOptions& sam_options);
};

} // namespace labw::art_modern