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

#ifndef ART_MODERN_LIBAM_SUPPORT_UTILS_HTS_UTILS_H
#define ART_MODERN_LIBAM_SUPPORT_UTILS_HTS_UTILS_H

#include "libam_support/utils/c_macro_utils.h"

#include <htslib/hts.h>

#include <stdlib.h>

ART_MODERN_BEGIN_C_DECLS

/**
 * @brief Know where we are.
 *
 * Note: no support for cram
 *
 * @param fp HTS file pointer
 * @return Where we are.
 */
size_t am_hts_tell(htsFile* fp);

ART_MODERN_END_C_DECLS

#endif // ART_MODERN_LIBAM_SUPPORT_UTILS_HTS_UTILS_H
