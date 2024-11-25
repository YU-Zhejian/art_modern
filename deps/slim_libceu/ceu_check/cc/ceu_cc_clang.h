#ifndef CEU_CC_CLANG_H
#define CEU_CC_CLANG_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif
#undef CEU_COMPILER_IS_CLANG

#if defined(__llvm__)
#define CEU_COMPILER_IS_LLVM
#endif

#if defined(__clang__)
#define CEU_COMPILER_IS_CLANG
#endif

#if defined(CEU_COMPILER_IS_CLANG)
#if !defined(CEU_COMPILER_NAME)
#define CEU_COMPILER_NAME "Clang"
#endif
#endif

#if defined(CEU_COMPILER_IS_LLVM)
#if !defined(CEU_COMPILER_NAME)
#define CEU_COMPILER_NAME "Unknown LLVM-based"
#endif
#endif
#endif // CEU_CC_CLANG_H
