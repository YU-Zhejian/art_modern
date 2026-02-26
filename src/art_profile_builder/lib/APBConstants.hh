/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include <cstdint>
namespace labw::art_modern {
enum class APB_FORMAT : std::uint8_t { FASTQ, SAM, BAM, CRAM, SRA, AUTO };
constexpr char APB_FORMAT_FASTQ_STR[] = "FASTQ";
constexpr char APB_FORMAT_SAM_STR[] = "SAM";
constexpr char APB_FORMAT_BAM_STR[] = "BAM";
constexpr char APB_FORMAT_CRAM_STR[] = "CRAM";
constexpr char APB_FORMAT_SRA_STR[] = "SRA";
constexpr char APB_FORMAT_AUTO_STR[] = "AUTO";
} // namespace labw::art_modern
