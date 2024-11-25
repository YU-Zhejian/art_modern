#ifndef CEU_OS_BSD_H
#define CEU_OS_BSD_H

#ifndef CEU_CHECK_OS_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_os_macro.h> instead!"
#endif
#if defined(__FreeBSD__)
#define CEU_ON_FreeBSD
#define CEU_PRIMARY_OS_TYPE "Free BSD"
#elif defined(__DragonFly__)
#define CEU_ON_DragonFlyBSD
#define CEU_PRIMARY_OS_TYPE "DragonFly BSD"
#elif defined(__NetBSD__)
#define CEU_ON_NetBSD
#define CEU_PRIMARY_OS_TYPE "NetBSD"
#elif defined(__OpenBSD__)
#define CEU_ON_OpenBSD
#define CEU_PRIMARY_OS_TYPE "OpenBSD"
#endif

#if defined(CEU_ON_OpenBSD) && defined(CEU_ON_FreeBSD) && defined(CEU_ON_DragonFlyBSD) && defined(CEU_ON_NetBSD)       \
    && defined(BSD) && defined(_BSD) && defined(__BSD) && defined(__BSD__) && defined(_BSD_)
#define CEU_ON_BSD
#ifndef CEU_PRIMARY_OS_TYPE
#define CEU_PRIMARY_OS_TYPE "Other BSD"
#endif
#endif

#endif
