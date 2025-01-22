/**
 * @brief Basic data types.
 */
#pragma once
#include <cstdint>

namespace labw::art_modern {
/** Phred quality score. */
using am_qual_t = uint8_t;
/** Cigar type. Should be 0b00 to 0b11. */
using am_cigar_type_t = uint8_t;
/** Cigar operation, length, or combined. */
using am_cigar_t = uint32_t;
/** ART quality distribution. */
using am_qual_dist_t = int;
/** Number of reads. */
using am_readnum_t = int64_t;
/** Length of files. */
using am_filelen_t = int64_t;
} // namespace labw::art_modern
