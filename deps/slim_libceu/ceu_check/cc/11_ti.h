#ifndef CEU_CC_TI_H
#define CEU_CC_TI_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif

/* Some TI compilers define __TI_COMPILER_VERSION without trailing underscores */
#if defined(__TI_COMPILER_VERSION__) || defined(__TI_COMPILER_VERSION)
#define CEU_COMPILER_IS_TI
#ifndef CEU_COMPILER_NAME
#if defined(__TI_ARM__) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__) /* TI ARMCL Compiler */
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_TI_ARMCL
#elif defined(__ARM_ARCH) && defined(__clang__) /* TI Arm Clang Compiler */
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_TI_ARMCLANG
#else
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_TI_UNKNOWN
#endif
#endif
#endif
#endif /* CEU_CC_TI_H */
