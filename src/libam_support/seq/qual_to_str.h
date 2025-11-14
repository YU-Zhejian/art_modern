//
// Created by yuzj on 11/1/25.
//

#ifndef ART_MODERN_LIBAM_SUPPORT_QUAL_TO_STR_H
#define ART_MODERN_LIBAM_SUPPORT_QUAL_TO_STR_H

#include "libam_support/Dtypes.h"
#include "libam_support/utils/c_macro_utils.h"

#include <stdlib.h>

ART_MODERN_BEGIN_C_DECLS

void qual_to_str_mmx(const am_qual_t* qual, char* str, size_t qlen);
void qual_to_str_sse2(const am_qual_t* qual, char* str, size_t qlen);
void qual_to_str_avx2(const am_qual_t* qual, char* str, size_t qlen);
void qual_to_str_for_loop(const am_qual_t* qual, char* str, size_t qlen);
void qual_to_str_comb(const am_qual_t* qual, char* str, size_t qlen);

ART_MODERN_END_C_DECLS

#endif // ART_MODERN_LIBAM_SUPPORT_QUAL_TO_STR_H
