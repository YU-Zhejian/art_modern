#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace labw::art_modern {
enum class CompressionType : std::uint8_t { NONE, GZIP, BGZIP, NOP };

static const std::vector<std::string> GZIP_EXTENSIONS = { ".gz", ".gzip", ".GZ" };
static const std::vector<std::string> VALID_COMPRESS_EXTENSIONS = { ".gz", ".bz2", ".xz", ".GZ", ".gzip", ".bgz" };
static const std::vector<std::string> BGZIP_EXTENSIONS = { ".bgz", ".BGZ" };

} // namespace labw::art_modern
