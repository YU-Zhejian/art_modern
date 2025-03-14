/**
 * @brief Basic data types.
 */
#pragma once
#include <cstdint>

namespace labw::art_modern {
/**
 * Phred quality score.
 *
 * This score will be involved in +/- operations so it must be signed.
 */
using am_qual_t = int8_t;
/** Cigar type. Should be 0b00 to 0b11. */
using am_cigar_type_t = uint8_t;
/**
 * Cigar operation, length, or combined.
 *
 * FIXME: This would cause error when doing +/- operations on cigar length
 * as they may create negative values.
 */
using am_cigar_t = uint32_t;
/** ART quality distribution. */
using am_qual_count_t = uint64_t;
/**
 * Number of reads.
 *
 * This MUST be a signed integer.
 * Otherwise -1 will be reverted to an extraordinarily large number.
 */
using am_readnum_t = int64_t;
/** Length of files. */
using am_filelen_t = int64_t;
} // namespace labw::art_modern
