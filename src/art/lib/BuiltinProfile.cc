/**
 * Copyright 2008-2016 Weichun Huang <whduke@gmail.com>
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "art/lib/BuiltinProfile.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <zlib.h>

#include <cstring>
#include <string>

namespace labw::art_modern {

namespace {
    constexpr int BUFF_SIZE = 4096;
    std::string decompress(const unsigned char* src, const std::size_t slen)
    {
        z_stream zs = {};

        // Initialize zlib stream for decompression
        if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
            BOOST_LOG_TRIVIAL(fatal) << "Failed to initialize zlib stream";
            abort_mpi();
        }

        zs.next_in = const_cast<unsigned char*>(src);
        zs.avail_in = slen;

        int ret = Z_OK;
        char outbuffer[BUFF_SIZE];
        std::string outstring;

        // Decompress the data
        do {
            zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
            zs.avail_out = sizeof(outbuffer);

            ret = inflate(&zs, 0);

            if (outstring.size() < zs.total_out) {
                outstring.append(outbuffer, zs.total_out - outstring.size());
            }
        } while (ret == Z_OK);

        // Clean up
        inflateEnd(&zs);

        // Check for errors
        if (ret != Z_STREAM_END) {
            BOOST_LOG_TRIVIAL(fatal) << "Error during zlib decompression: " << zs.msg;
            abort_mpi();
        }

        return outstring;
    }
} // namespace

BuiltinProfile::BuiltinProfile(const unsigned char* profile1_encoded, const std::size_t profile1_encoded_size,
    const unsigned char* profile2_encoded, const std::size_t profile2_encoded_size)
    : r1_profile(decompress(profile1_encoded, profile1_encoded_size))
    , r2_profile(profile2_encoded[0] != '\0' ? decompress(profile2_encoded, profile2_encoded_size) : "")
{
}
} // namespace labw::art_modern
