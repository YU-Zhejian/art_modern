/*!
 * @file ceu_check_cxx_std_macro.hh
 * @brief Check Compile-Time C++ standard.
 * This function defines macros only.
 */

#pragma once

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
#elif CEU_CXX_STD_VERSION_MACRO == 199711L
#define CEU_CXX_STD "98"
#elif CEU_CXX_STD_VERSION_MACRO == 201103L
#define CEU_CXX_STD "11"
#elif CEU_CXX_STD_VERSION_MACRO == 201402L
#define CEU_CXX_STD "14"
#elif CEU_CXX_STD_VERSION_MACRO == 201703L
#define CEU_CXX_STD "17"
#elif CEU_CXX_STD_VERSION_MACRO == 202002L
#define CEU_CXX_STD "20"
/* GCC 11 to 13 report 202100L for -std=c++23 and -std=c++2b; 202302L is the ratified value. */
#elif CEU_CXX_STD_VERSION_MACRO == 202100L || CEU_CXX_STD_VERSION_MACRO == 202302L
#define CEU_CXX_STD "23"
/*
 * C++26 is not ratified yet, so no final value is assigned. GCC 14+ and Clang 17+ report the placeholder 202400L
 * for -std=c++26 and -std=c++2c. Every ratified revision has used YYYYMM of its publication, so accept any 2026
 * month too, in order to keep reporting "26" once compilers switch to the final value.
 */
#elif CEU_CXX_STD_VERSION_MACRO == 202400L                                                                             \
    || (CEU_CXX_STD_VERSION_MACRO >= 202601L && CEU_CXX_STD_VERSION_MACRO <= 202612L)
#define CEU_CXX_STD "26"
#elif CEU_CXX_STD_VERSION_MACRO > 202612L
#define CEU_CXX_STD "post-26"
#else
#define CEU_CXX_STD "unknown"
#endif
