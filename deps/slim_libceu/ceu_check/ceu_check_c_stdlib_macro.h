/**
 * TODO: Check libc
 */
#ifndef CEU_CHECK_C_STDLIB_MACRO_H
#define CEU_CHECK_C_STDLIB_MACRO_H

#include "libceu_stddef.h"

#if defined(CEU_HAVE_INCLUDE_FEATURES_H) && (CEU_HAVE_INCLUDE_FEATURES_H == 1)
#include <features.h>
#endif

#include "ceu_check/ceu_check_cc_macro.h"

#include <limits.h> /* NOLINT for __GLIBC__ */

#if defined(__ANDROID__)
#define CEU_LIBC_IS_ANDROID
#elif defined (CEU_COMPILER_IS_MSVC)
/**
 * Should be either MSVCRT or UCRT.
 */
#define CEU_LIBC_IS_MSVC
#else
#if defined(__UCLIBC__)
#define CEU_LIBC_IS_UCLIBC
#elif defined(__dietlibc__)
#define CEU_LIBC_IS_DIETLIBC
#elif defined(__GLIBC__) || defined(__GNU_LIBRARY__)
#define CEU_LIBC_IS_GNU
#else
#include <stdarg.h>
/* This part of code is from GNU config.guess */
/* First heuristic to detect musl libc.  */
#ifdef __DEFINED_va_list
#define CEU_LIBC_IS_MUSL
#else
#define CEU_LIBC_IS_UNKNOWN
#endif
#endif
#endif


#endif /* CEU_CHECK_C_STDLIB_MACRO_H */
