/**
 * TODO: Check libc
 */

#include <limits.h> /* NOLINT for __GLIBC__ */

#if defined(__ANDROID__)
#define CEU_LIBC_IS_ANDROID
#else
#include <features.h>
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
