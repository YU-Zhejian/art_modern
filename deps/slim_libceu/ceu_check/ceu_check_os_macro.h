/**
 * @file ceu_check_os_macro.h
 * @author YU Zhejian
 * @brief Check operating system pre-defined macros.
 * @version 0.1
 * @date 2024-05-13
 *
 */
#ifndef CEU_CHECK_OS_MACRO_H
#define CEU_CHECK_OS_MACRO_H

#ifdef CEU_PRIMARY_OS_TYPE
#undef CEU_PRIMARY_OS_TYPE
#endif

/* Inclusion order of following lines should be preserved to avoid conflicts. */
/* clang-format off */
#include <ceu_check/os/ceu_os_mainframe_unix.h> /* NOLINT */
#include <ceu_check/os/ceu_os_bsd.h> /* NOLINT */
#include <ceu_check/os/ceu_os_windows.h> /* NOLINT */
/* clang-format on */
#if defined(__HAIKU__)
#define CEU_ON_HAIKU
#define CEU_PRIMARY_OS_TYPE "Haiku"
#endif /* __HAIKU__ */

#if defined(__ANDROID__)
#define CEU_ON_ANDROID
#define CEU_PRIMARY_OS_TYPE "Android"
#endif /* __ANDROID__ */

#if defined(__MACH__)
#define CEU_ON_MACH
#if defined(__APPLE__)
#define CEU_ON_APPLE
#define CEU_PRIMARY_OS_TYPE "Apple MACH"
#elif defined(__gnu_hurd__)
#define CEU_ON_GNU_HURD
#define CEU_PRIMARY_OS_TYPE "GNU Hurd"
#else
#define CEU_PRIMARY_OS_TYPE "MACH UNKNOWN"
#endif
#endif /* __MACH__ */

#if defined(__linux) || defined(__linux__) || defined(linux) || defined(__gnu_linux__)
#define CEU_ON_GNU_LINUX
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "GNU/Linux" /* NOLINT */
#endif
#endif /* __linux */

#if defined(__unix__) || defined(__unix) || defined(unix) || defined(_unix) || defined(_unix_) || defined(CEU_ON_MACH) \
    || defined(CEU_ON_GNU_LINUX) || defined(CEU_ON_BSD) || defined(CEU_ON_SOLARIS) || defined(CEU_ON_HP_UX)            \
    || defined(CEU_ON_AIX)
#define CEU_ON_UNIX
#endif /* __unix__ */

/* Note that MSYS and MinGW are NOT POSIX. */
#if defined(CEU_ON_UNIX) || defined(CEU_ON_CYGWIN)
#define CEU_ON_POSIX
#endif

#if defined(__WINE__)
#define CEU_ON_WINE
#endif

#ifndef CEU_PRIMARY_OS_TYPE
#ifdef CEU_ON_UNIX
#define CEU_PRIMARY_OS_TYPE "Other UNIX"
#else
#define CEU_PRIMARY_OS_TYPE "UNKNOWN"
#endif
#endif /* CEU_PRIMARY_OS_TYPE */

#endif /* CEU_CHECK_OS_MACRO_H */
