#ifndef CEU_CC_AMD_CLANG_H
#define CEU_CC_AMD_CLANG_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif
#if defined(__aocc__)
#define CEU_COMPILER_IS_AOCC
#endif
#if (defined(CEU_COMPILER_IS_AOCC) && !defined(CEU_COMPILER_NAME))
#define CEU_COMPILER_NAME "AMD Optimizing C/C++ and Fortran Compilers (AOCC)"
#endif
#endif
