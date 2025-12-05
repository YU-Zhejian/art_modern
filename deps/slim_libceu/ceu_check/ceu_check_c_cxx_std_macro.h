/*!
 * @file ceu_check_c_cxx_std_macro.h
 * @brief Check Compile-Time C/C++ standard.
 * This function defines macros only.
 */

#ifndef CEU_CHECK_C_CXX_STD_MACRO_H
#define CEU_CHECK_C_CXX_STD_MACRO_H

#include "ceu_check/ceu_constants.h"

#ifdef CEU_CXX_STD_VERSION_MACRO
#undef CEU_CXX_STD_VERSION_MACRO
#endif
#if defined(_MSVC_LANG)
#define CEU_CXX_STD_VERSION_MACRO _MSVC_LANG
#elif defined(__cplusplus)
#define CEU_CXX_STD_VERSION_MACRO __cplusplus
#endif

#ifndef CEU_CXX_STD_VERSION_MACRO
#define CEU_CXX_STD CEU_UNDEFINED
#elif CEU_CXX_STD_VERSION_MACRO < 199711L
#define CEU_CXX_STD "pre-98"
#elif CEU_CXX_STD_VERSION_MACRO >= 202100L
#define CEU_CXX_STD "post-23"
#elif CEU_CXX_STD_VERSION_MACRO == 202100L
#define CEU_CXX_STD "23"
#elif CEU_CXX_STD_VERSION_MACRO == 202002L
#define CEU_CXX_STD "20"
#elif CEU_CXX_STD_VERSION_MACRO == 201703L
#define CEU_CXX_STD "17"
#elif CEU_CXX_STD_VERSION_MACRO == 201402L
#define CEU_CXX_STD "14"
#elif CEU_CXX_STD_VERSION_MACRO == 201103L
#define CEU_CXX_STD "11"
#elif CEU_CXX_STD_VERSION_MACRO == 199711L
#define CEU_CXX_STD "98"
#endif

#ifdef CEU_C_STD_VERSION_MACRO
#undef CEU_C_STD_VERSION_MACRO
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

#endif /* CEU_CHECK_C_CXX_STD_MACRO_H  */
