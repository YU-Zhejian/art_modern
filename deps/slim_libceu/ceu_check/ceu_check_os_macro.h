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

#ifdef CEU_UNDER_DOXYGEN

/*!
 *@def CEU_ON_MSYS
 *@brief Indicating CEU is compiled on MSYS
 */
#define CEU_ON_MSYS

/*!
 *@def CEU_ON_CYGWIN
 *@brief Indicating CEU is compiled on CYGWin.
 */
#define CEU_ON_CYGWIN

/*!
 *@def CEU_ON_MINGW64
 *@brief Indicating CEU is compiled on MinGW 64.
 */
#define CEU_ON_MINGW64

/*!
 *@def CEU_ON_MINGW32
 *@brief Indicating CEU is compiled on MinGW 32.
 */
#define CEU_ON_MINGW32

/*!
 *@def CEU_ON_CYGWIN_LIKE
 *@brief Indicating CEU is compiled on CygWin/MSYS/MinGW or friends.
 */
#define CEU_ON_CYGWIN_LIKE

/*!
 *@def CEU_ON_WINDOWS_64
 *@brief Indicating CEU is compiled on 64-bit Microsoft Windows.
 */
#define CEU_ON_WINDOWS_64

/*!
 *@def CEU_ON_WINDOWS_32
 *@brief Indicating CEU is compiled on 32-bit Microsoft Windows.
 */
#define CEU_ON_WINDOWS_32

/*!
 *@def CEU_ON_WINDOWS
 *@brief Indicating CEU is compiled on Microsoft Windows.
 */
#define CEU_ON_WINDOWS

/*!
 *@def CEU_ON_FreeBSD
 *@brief Indicating CEU is compiled on Free BSD.
 */
#define CEU_ON_FreeBSD

/*!
 *@def CEU_ON_DragonFlyBSD
 *@brief Indicating CEU is compiled on DragonFly BSD.
 */
#define CEU_ON_DragonFlyBSD

/*!
 *@def CEU_ON_NetBSD
 *@brief Indicating CEU is compiled on NetBSD.
 */
#define CEU_ON_NetBSD

/*!
 *@def CEU_ON_OpenBSD
 *@brief Indicating CEU is compiled on OpenBSD.
 */
#define CEU_ON_OpenBSD

/*!
 *@def CEU_ON_BSD
 *@brief Indicating CEU is compiled on some BSD.
 */
#define CEU_ON_BSD

/*!
 *@def CEU_ON_SOLARIS
 *@brief Indicating CEU is compiled on Solaris.
 */
#define CEU_ON_SOLARIS

/*!
 *@def CEU_ON_HP_UX
 *@brief Indicating CEU is compiled on HP-UX.
 */
#define CEU_ON_HP_UX

/*!
 *@def CEU_ON_AIX
 *@brief Indicating CEU is compiled on IBM AIX.
 */
#define CEU_ON_AIX

/*!
 *@def CEU_ON_HAIKU
 *@brief Indicating CEU is compiled on Haiku.
 */
#define CEU_ON_HAIKU

/*!
 *@def CEU_ON_ANDROID
 *@brief Indicating CEU is compiled on Android.
 */
#define CEU_ON_ANDROID

/*!
 *@def CEU_PRIMARY_OS_TYPE
 *@brief A human-readable operating system name.
 */
#define CEU_PRIMARY_OS_TYPE /* implementation-defined */

/*!
 *@def CEU_ON_MACH
 *@brief Indicating CEU is compiled on MACH macro kernel.
 */
#define CEU_ON_MACH

/*!
 *@def CEU_ON_APPLE
 *@brief Indicating CEU is compiled on Apple.
 */
#define CEU_ON_APPLE

/*!
 *@def CEU_ON_GNU_HURD
 *@brief Indicating CEU is compiled on GNU Hurd.
 */
#define CEU_ON_GNU_HURD

/*!
 *@def CEU_ON_GNU_LINUX
 *@brief Indicating CEU is compiled on GNU/Linux.
 */
#define CEU_ON_GNU_LINUX

#define CEU_ON_UNIX

#define CEU_ON_POSIX

#define CEU_ON_WINE

#else

#ifdef CEU_PRIMARY_OS_TYPE
#undef CEU_PRIMARY_OS_TYPE
#endif

// Inclusion order of following lines should be preserved
// clang-format off
#include <ceu_check/os/ceu_os_mainframe_unix.h> // NOLINT
#include <ceu_check/os/ceu_os_bsd.h> // NOLINT
#include <ceu_check/os/ceu_os_windows.h> // NOLINT

// clang-format on
#if defined(__HAIKU__)
#define CEU_ON_HAIKU
#define CEU_PRIMARY_OS_TYPE "Haiku"
#endif

#if defined(__ANDROID__)
#define CEU_ON_ANDROID
#define CEU_PRIMARY_OS_TYPE "Android"
#endif

#if defined(__MACH__)
#define CEU_ON_MACH
#if defined(__APPLE__)
#define CEU_ON_APPLE
#elif defined(__gnu_hurd__)
#define CEU_ON_GNU_HURD
#define CEU_PRIMARY_OS_TYPE "GNU Hurd"
#else
#define CEU_PRIMARY_OS_TYPE "MACH UNKNOWN"
#endif
#endif

#if defined(__linux) || defined(__linux__) || defined(linux) || defined(__gnu_linux__)
#define CEU_ON_GNU_LINUX
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "GNU/Linux"
#endif
#endif

#if defined(__unix__) || defined(__unix) || defined(unix) || defined(_unix) || defined(_unix_) || defined(CEU_ON_MACH) \
    || defined(CEU_ON_GNU_LINUX) || defined(CEU_ON_BSD) || defined(CEU_ON_SOLARIS) || defined(CEU_ON_HP_UX)            \
    || defined(CEU_ON_AIX)
#define CEU_ON_UNIX
#endif


// Note that MSYS and MinGW are NOT POSIX.
#if defined(CEU_ON_UNIX) || defined(CEU_ON_CYGWIN_LIKE)
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
#endif

#endif
#endif // CEU_CHECK_OS_MACRO_H
