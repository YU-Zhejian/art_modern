/**
 * TODO: From <https://www.boost.org/doc/libs/1_64_0/boost/functional/hash/hash.hpp>
 */

#pragma once

namespace labw::art_modern {
template <typename SizeT> void boost_hash_combine_impl(SizeT& seed, SizeT value)
{
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

} // namespace labw::art_modern
