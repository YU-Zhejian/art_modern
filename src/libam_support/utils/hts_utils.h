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

#include "libam_support/Dtypes.h"
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

/*! @function
 @abstract  Get the CIGAR array
 @param  b  pointer to an alignment
 @return    pointer to the CIGAR array

 @discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
 */
#ifdef __cplusplus
#define am_bam_get_cigar(b) (reinterpret_cast<am_cigar_t*>((b)->data + (b)->core.l_qname))
#else
#define am_bam_get_cigar(b) ((am_cigar_t*)((b)->data + (b)->core.l_qname))
#endif

ART_MODERN_END_C_DECLS

#endif /** ART_MODERN_LIBAM_SUPPORT_UTILS_HTS_UTILS_H */
