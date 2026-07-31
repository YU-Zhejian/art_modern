/**
 * Copyright 2024-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
