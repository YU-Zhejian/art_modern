/* Include this file before "ceu_cc_clang.h"!   */

#ifndef CEU_CC_ICC_H
#define CEU_CC_ICC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ECC)
#define CEU_COMPILER_IS_ICC
#if defined(__llvm__)
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_ICC_LLVM
#else
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_ICC
#endif
#endif
#endif /* CEU_CC_ICC_H */
