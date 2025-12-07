/*
 * TODO: Check TI armcl vs. tiarmclang
 */

#ifndef CEU_CC_TI_H
#define CEU_CC_TI_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif

/* Some TI compilers define __TI_COMPILER_VERSION without trailing underscores */
#if defined(__TI_COMPILER_VERSION__) || defined(__TI_COMPILER_VERSION)
#define CEU_COMPILER_IS_TI
#ifndef CEU_COMPILER_NAME
#define CEU_COMPILER_NAME CEU_COMPILER_NMAME_TI
#endif
#endif
#endif /* CEU_CC_TI_H */
