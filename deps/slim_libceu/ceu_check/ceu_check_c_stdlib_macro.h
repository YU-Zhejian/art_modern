#ifndef CEU_CHECK_C_STDLIB_MACRO_H
#define CEU_CHECK_C_STDLIB_MACRO_H

#include "libceu_stddef.h"

#if defined(CEU_HAVE_INCLUDE_FEATURES_H) && (CEU_HAVE_INCLUDE_FEATURES_H == 1)
#include <features.h> /* NOLINT*/
#endif

#if defined(CEU_HAVE_INCLUDE_SYS_CDEFS_H) && (CEU_HAVE_INCLUDE_SYS_CDEFS_H == 1)
#include <sys/cdefs.h> /* NOLINT*/
#endif

#include "ceu_check/ceu_check_cc_macro.h"

#include <limits.h> /* NOLINT for __GLIBC__ */

#undef CEU_LIBC_NAME

#if defined(__ANDROID__)
#define CEU_LIBC_IS_ANDROID
#define CEU_LIBC_NAME "Bionic"
#elif defined(CEU_COMPILER_IS_MSVC)
/**
 * Should be either MSVCRT or UCRT.
 */
#define CEU_LIBC_IS_MSVC
#define CEU_LIBC_NAME "MSVCRT or UCRT"
#elif defined(__UCLIBC__)
#define CEU_LIBC_IS_UCLIBC
#define CEU_LIBC_NAME "uClibc"
#elif defined(__dietlibc__)
#define CEU_LIBC_IS_DIETLIBC
#define CEU_LIBC_NAME "dietlibc"
#elif defined(__NEWLIB__)
/* Newlib is also used in Cygwin */
#define CEU_LIBC_IS_NEWLIB
#define CEU_LIBC_NAME "Newlib"
#elif defined(__GLIBC__) || defined(__GNU_LIBRARY__)
#define CEU_LIBC_IS_GNU
#define CEU_LIBC_NAME "GNU C Library (glibc)"
#else
#include <stdarg.h>
/* This part of the code is from GNU config.guess */
/* First heuristic to detect musl libc.  */
#ifdef __DEFINED_va_list
#define CEU_LIBC_IS_MUSL
#define CEU_LIBC_NAME "(Possibly) musl libc"
#else
#define CEU_LIBC_IS_UNKNOWN
#define CEU_LIBC_NAME "Unknown C Standard Library"
#endif
#endif

#endif /* CEU_CHECK_C_STDLIB_MACRO_H */
