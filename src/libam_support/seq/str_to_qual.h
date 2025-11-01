//
// Created by yuzj on 11/1/25.
//

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

ART_MODERN_END_C_DECLS

#endif //ART_MODERN_LIBAM_SUPPORT_SEQ_STR_TO_QUAL_H
