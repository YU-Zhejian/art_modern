#ifndef CEU_CC_INTEL_CLANG_H
#define CEU_CC_INTEL_CLANG_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif
#if defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define CEU_COMPILER_IS_INTEL_CLANG
#define CEU_COMPILER_NAME "Intel(R) oneAPI DPC++/C++ Compiler (ICX)"
#endif
#endif
