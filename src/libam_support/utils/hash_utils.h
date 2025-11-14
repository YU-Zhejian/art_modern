/**
 *
 * From <https://www.boost.org/doc/libs/1_64_0/boost/functional/hash/hash.hpp>
 *
 * Original license:
 *
 * Copyright 2005-2014 Daniel James.
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Based on Peter Dimov's proposal
 *  http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
 *  issue 6.18.
 *
 *  This also contains public domain code from MurmurHash. From the
 *  MurmurHash header:
 *
 * MurmurHash3 was written by Austin Appleby, and is placed in the public
 * domain. The author hereby disclaims copyright to this source code.
 */

#ifndef LABW_ART_MODERN_LIBAM_SUPPORT_UTILS_HASH_UTILS_H
#define LABW_ART_MODERN_LIBAM_SUPPORT_UTILS_HASH_UTILS_H
#include "libam_support/utils/c_macro_utils.h"

#include <stdint.h>

ART_MODERN_BEGIN_C_DECLS

inline void am_boost_hash_combine_impl32(uint32_t& seed, uint32_t value)
{
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Improved by Tongyi Lingma
inline void am_boost_hash_combine_impl64(uint64_t& seed, uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 13) + (seed >> 19);
}

ART_MODERN_END_C_DECLS

#endif // LABW_ART_MODERN_LIBAM_SUPPORT_UTILS_HASH_UTILS_H
