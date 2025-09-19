/**
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

#pragma once

#include <htslib/hts.h>

#include <cstddef>

namespace labw::art_modern {
/**
 * @brief Know where we are.
 *
 * Note: no support for cram
 *
 * @param fp HTS file pointer
 * @return Where we are.
 */
std::size_t hts_tell(htsFile* fp);
} // namespace labw::art_modern
