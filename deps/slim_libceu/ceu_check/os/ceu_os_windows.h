#ifndef CEU_OS_WINDOWS_H
#define CEU_OS_WINDOWS_H

#ifndef CEU_CHECK_OS_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_os_macro.h> instead!"
#endif
#if defined(__MSYS__)
#define CEU_PRIMARY_OS_TYPE "MSYS"

#define CEU_ON_MSYS
#endif

#if defined(__CYGWIN__)
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "CYGWIN"
#endif
#define CEU_ON_CYGWIN
#endif

#if defined(__MINGW64__)
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "MINGW64"
#endif
#define CEU_ON_MINGW64
#endif

#if defined(__MINGW32__) || defined(__WINNT__) || defined(__WINNT)
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "MINGW32"
#endif
#define CEU_ON_MINGW32
#endif

#if defined(CEU_ON_CYGWIN) || defined(CEU_ON_MINGW64) || defined(CEU_ON_MINGW32) || defined(CEU_ON_MSYS)
#define CEU_ON_CYGWIN_LIKE
#endif

#if defined(WIN64) || defined(_WIN64) || defined(__WIN64) || defined(__WIN64__)
#define CEU_ON_WINDOWS_64
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "Microsoft Windows 64"
#endif
#elif (defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(__WIN32__))
#define CEU_ON_WINDOWS_32
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "Microsoft Windows 32"
#endif
#endif

#if defined(CEU_ON_WINDOWS_32) || defined(CEU_ON_WINDOWS_64)
#define CEU_ON_WINDOWS
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "Microsoft Windows"
#endif
#endif

#endif
