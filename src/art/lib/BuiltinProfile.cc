#include "art/lib/BuiltinProfile.hh"

// #include <base64.h>

#include <zlib.h>

#include <cstring>
#include <stdexcept>
#include <string>

namespace labw::art_modern {

namespace {
    constexpr int BUFF_SIZE = 4096;
    std::string decompress(const unsigned char* src, std::size_t slen)
    {
        z_stream zs;
        std::memset(&zs, 0, sizeof(zs));

        // Initialize zlib stream for decompression
        if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
            throw std::runtime_error("Failed to initialize zlib stream");
        }

        zs.next_in = const_cast<unsigned char*>(src);
        zs.avail_in = slen;

        int ret;
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
            throw std::runtime_error("Error during zlib decompression: " + std::string(zs.msg));
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
