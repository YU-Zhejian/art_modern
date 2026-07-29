#pragma once

#include <array>
#include <cstdint>
#include <string_view>

namespace labw::art_modern {
enum class CompressionType : std::uint8_t { NONE, GZIP, BGZIP, NOP };

// Declared constexpr so that they carry no dynamic initialization. A namespace-scope `static const
// std::vector<std::string>` in a header has internal linkage, which means one heap-allocating copy per
// translation unit that includes this header, constructed before main().
constexpr std::array<std::string_view, 3> GZIP_EXTENSIONS { ".gz", ".gzip", ".GZ" };
constexpr std::array<std::string_view, 6> VALID_COMPRESS_EXTENSIONS { ".gz", ".bz2", ".xz", ".GZ", ".gzip", ".bgz" };
constexpr std::array<std::string_view, 2> BGZIP_EXTENSIONS { ".bgz", ".BGZ" };

} // namespace labw::art_modern
