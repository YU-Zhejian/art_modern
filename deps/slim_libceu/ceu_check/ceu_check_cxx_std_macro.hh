/*!
 * @file ceu_check_cxx_std_macro.hh
 * @brief Check Compile-Time C++ standard.
 * This function defines macros only.
 */

#pragma once

#include "ceu_check/ceu_constants.h" /* NOLINT: for CEU_UNDEFINED */

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
