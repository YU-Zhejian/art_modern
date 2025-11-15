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

#ifndef ART_MODERN_LIBAM_SUPPORT_SEQ_STR_TO_QUAL_H
#define ART_MODERN_LIBAM_SUPPORT_SEQ_STR_TO_QUAL_H

#include "libam_support/Dtypes.h"
#include "libam_support/utils/c_macro_utils.h"

#include <stdlib.h>

ART_MODERN_BEGIN_C_DECLS

void str_to_qual_avx2(am_qual_t* qual, const char* str, size_t qlen);
void str_to_qual_mmx(am_qual_t* qual, const char* str, size_t qlen);
void str_to_qual_sse2(am_qual_t* qual, const char* str, size_t qlen);
void str_to_qual_for_loop(am_qual_t* qual, const char* str, size_t qlen);
void str_to_qual_comb(am_qual_t* qual, const char* str, size_t qlen);

ART_MODERN_END_C_DECLS

#endif // ART_MODERN_LIBAM_SUPPORT_SEQ_STR_TO_QUAL_H
