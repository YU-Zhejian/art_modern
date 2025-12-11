/*!
 * @file ceu_check_c_std_macro.h
 * @brief Check Compile-Time C standard.
 * This function defines macros only.
 */

#ifndef CEU_CHECK_C_STD_MACRO_H
#define CEU_CHECK_C_STD_MACRO_H

#include "ceu_utils/ceu_constants.h" /* NOLINT: for CEU_UNDEFINED */

#if defined(__cplusplus)
#error "This is a C-only header. Do not include it in C++ code!"
#endif

#if defined(__STDC_VERSION__)
#define CEU_C_STD_VERSION_MACRO __STDC_VERSION__
#endif

#ifndef CEU_C_STD_VERSION_MACRO
#ifdef __STDC__
#define CEU_C_STD "unknown ISO C"
#else
#define CEU_C_STD CEU_UNDEFINED
#endif
#elif CEU_C_STD_VERSION_MACRO < 199409L
#define CEU_C_STD "pre-94"
#elif CEU_C_STD_VERSION_MACRO == 199409L
#define CEU_C_STD "94"
#elif CEU_C_STD_VERSION_MACRO == 199901L
#define CEU_C_STD "99"
#elif CEU_C_STD_VERSION_MACRO == 201112L
#define CEU_C_STD "11"
#elif CEU_C_STD_VERSION_MACRO == 201710L
#define CEU_C_STD "17"
#elif CEU_C_STD_VERSION_MACRO == 202311L
#define CEU_C_STD "23"
#elif CEU_C_STD_VERSION_MACRO > 202311L
#define CEU_C_STD "post-23"
#else
#define CEU_C_STD "unknown"
#endif

#endif /* CEU_CHECK_C_STD_MACRO_H  */
