#ifndef CEU_CC_EDG_H
#define CEU_CC_EDG_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif

#if defined(__EDG__)
#define CEU_COMPILER_IS_EDG
// Do not define compiler name since there are multiple compilers that uses EDG frontend.
#endif
#endif
