/**
 * Copyright 2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

#ifndef ART_MODERN_LIBAM_SUPPORT_UTILS_C_MACROS_H
#define ART_MODERN_LIBAM_SUPPORT_UTILS_C_MACROS_H

// See: https://www.gnu.org/software/libtool/manual/html_node/C-header-files.html
#undef ART_MODERN_BEGIN_C_DECLS
#undef ART_MODERN_END_C_DECLS

#ifdef __cplusplus
#define ART_MODERN_BEGIN_C_DECLS extern "C" {
#define ART_MODERN_END_C_DECLS }
#else
#define ART_MODERN_BEGIN_C_DECLS /* empty */
#define ART_MODERN_END_C_DECLS /* empty */
#endif
#if defined(__SSE2__) || defined(__AVX2__) || defined(__MMX__)
#ifdef _MSC_VER
#define INCLUDE_INTEL_SIMD <intrin.h>
#else
#define INCLUDE_INTEL_SIMD <x86intrin.h>
#endif
#else
#undef INCLUDE_INTEL_SIMD
#endif

#endif // ART_MODERN_LIBAM_SUPPORT_UTILS_C_MACROS_H
