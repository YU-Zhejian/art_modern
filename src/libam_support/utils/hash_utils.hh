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

#pragma once

namespace labw::art_modern {
template <typename SizeT> void boost_hash_combine_impl(SizeT& seed, SizeT value)
{
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

} // namespace labw::art_modern
