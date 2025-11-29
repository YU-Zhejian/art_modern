#ifndef CEU_CC_GCC_H
#define CEU_CC_GCC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif

#if defined(__GNUC__)
#define CEU_COMPILER_IS_GCC
#if !defined(CEU_COMPILER_NAME)
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_GCC
#endif
#endif
#endif /* CEU_CC_GCC_H */
