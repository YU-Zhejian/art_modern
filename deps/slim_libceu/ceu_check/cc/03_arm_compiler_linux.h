#ifndef CEU_CC_ARM_COMPILER_LINUX_H
#define CEU_CC_ARM_COMPILER_LINUX_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif
#if defined(__ARM_LINUX_COMPILER__)
#define CEU_COMPILER_IS_ARM_COMPILER_LINUX
#endif
#if (defined(CEU_COMPILER_IS_ARM_COMPILER_LINUX) && !defined(CEU_COMPILER_NAME))
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_ARM_COMPILER_LINUX
#endif
#endif /* CEU_CC_ARM_COMPILER_LINUX_H */
