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
#ifndef LIBAM_SUPPORT_CONSTANTS_H
#define LIBAM_SUPPORT_CONSTANTS_H

// This is a C header; the macro is consumed from C sources (see src/libam_support/seq/str_to_qual.c, where it is
// cast and fed to SIMD intrinsics), so it cannot become a C++ enum.
// NOLINTNEXTLINE(modernize-macro-to-enum,cppcoreguidelines-macro-to-enum,cppcoreguidelines-macro-usage)
#define AM_PHRED_OFFSET 33

#endif // LIBAM_SUPPORT_CONSTANTS_H
