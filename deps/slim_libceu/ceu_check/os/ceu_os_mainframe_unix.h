#ifndef CEU_OS_MAINFRAME_UNIX_H
#define CEU_OS_MAINFRAME_UNIX_H

#ifndef CEU_CHECK_OS_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_os_macro.h> instead!"
#endif

#if defined(__sun) && defined(__SVR4)
#define CEU_ON_SOLARIS
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "Solaris"
#endif
#endif

#if defined(__hpux)
#define CEU_ON_HP_UX
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "HP-UX"
#endif
#endif

#if defined(_AIX)
#define CEU_ON_AIX
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "IBM AIX"
#endif
#endif

#endif
