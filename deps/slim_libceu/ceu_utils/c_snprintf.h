#ifndef CEU_CHECK_C_SNPRINTF_H
#define CEU_CHECK_C_SNPRINTF_H

#include "ceu_check/ceu_check_c_cxx_compiler_macro.h" /* NOLINT: For CEU_COMPILER_IS_MSVC */

#include <stdio.h> /* NOLINT: For snprintf or _snprintf_s */

#ifdef CEU_COMPILER_IS_MSVC
/** Macro for _snprintf_s with truncation behavior similar to snprintf  */
#define CEU_SNPRINTF(buf, buf_size, ...) _snprintf_s(buf, buf_size, _TRUNCATE, __VA_ARGS__)
#else
/** Macro for snprintf */
#define CEU_SNPRINTF(buf, buf_size, ...) snprintf(buf, buf_size, __VA_ARGS__)
#endif

#endif /* CEU_CHECK_C_SNPRINTF_H */
