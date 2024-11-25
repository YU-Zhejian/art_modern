#ifndef CEU_CC_TURBOC_H
#define CEU_CC_TURBOC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif

#if defined(__BORLANDC__) || defined(__CODEGEARC__) || defined(__BCPLUSPLUS__) || defined(__TCPLUSPLUS__)              \
    || defined(__TURBOC__)
#define CEU_COMPILER_IS_BORLAND
#if defined(__clang__)
#define CEU_COMPILER_NAME "Borland C++Builder Compilers with Clang support"
#else
#define CEU_COMPILER_NAME "Borland C++Builder Compilers"
#endif

#undef CEU_CPPB_VERSION
#if __TURBOC__ < 0x0570
#define CEU_CPPB_VERSION "unknown before Developer Studio 2006"
#elif __TURBOC__ == 0x0570
#define CEU_CPPB_VERSION "Developer Studio 2006"
#elif __TURBOC__ == 0x0590
#define CEU_CPPB_VERSION "C++Builder 2007"
#elif __TURBOC__ == 0x0591
#define CEU_CPPB_VERSION "update 1 to C++Builder 2007"
#elif __TURBOC__ == 0x0592
#define CEU_CPPB_VERSION "RAD Studio 2007"
#elif __TURBOC__ == 0x0593
#define CEU_CPPB_VERSION "December update to RAD Studio 2007"
#elif __TURBOC__ == 0x0610
#define CEU_CPPB_VERSION "C++Builder 2009 (with/without Update 1)"
#elif __TURBOC__ == 0x0620
#define CEU_CPPB_VERSION "C++Builder 2010 (with/without Update 1)"
#elif __TURBOC__ == 0x0621
#define CEU_CPPB_VERSION "C++Builder 2010 Update 2"
#elif __TURBOC__ == 0x0630
#define CEU_CPPB_VERSION "C++Builder XE"
#elif __TURBOC__ == 0x0631
#define CEU_CPPB_VERSION "C++Builder XE Update 1"
#elif __TURBOC__ == 0x0640
#define CEU_CPPB_VERSION "C++Builder XE2"
#elif __TURBOC__ == 0x0650
#define CEU_CPPB_VERSION "C++Builder XE3 (with/without Update 1)"
#elif __TURBOC__ == 0x0651
#define CEU_CPPB_VERSION "C++Builder XE3 January 2013 update"
#elif __TURBOC__ == 0x0660
#define CEU_CPPB_VERSION "C++Builder XE4"
#elif __TURBOC__ == 0x0670
#define CEU_CPPB_VERSION "C++Builder XE5"
#elif __TURBOC__ == 0x0680
#define CEU_CPPB_VERSION "C++Builder XE6"
#elif __TURBOC__ == 0x0690
#define CEU_CPPB_VERSION "C++Builder XE7"
#elif __TURBOC__ == 0x0700
#define CEU_CPPB_VERSION "C++Builder XE8"
#elif __TURBOC__ == 0x0710
#define CEU_CPPB_VERSION "C++Builder Seattle"
#elif __TURBOC__ == 0x0711
#define CEU_CPPB_VERSION "C++Builder Seattle Subscription Update 1"
#elif __TURBOC__ == 0x0720
#define CEU_CPPB_VERSION "C++Builder Berlin"
#elif __TURBOC__ == 0x0730
#define CEU_CPPB_VERSION "C++Builder Tokyo"
#elif __TURBOC__ == 0x0740
#define CEU_CPPB_VERSION "C++Builder Rio"
#elif __TURBOC__ == 0x0750
#define CEU_CPPB_VERSION "C++Builder Sydney"
#elif __TURBOC__ == 0x0760
#define CEU_CPPB_VERSION "C++Builder Alexandria"
#elif __TURBOC__ > 0x0760
#define CEU_CPPB_VERSION "unknown after C++Builder Alexandria"
#endif
#endif
#endif
