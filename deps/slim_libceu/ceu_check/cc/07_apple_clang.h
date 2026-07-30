#ifndef CEU_CC_APPLE_CLANG_H
#define CEU_CC_APPLE_CLANG_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif

/*
 * Apple Clang is the Clang shipped with Xcode and the Apple Command Line Tools. It defines
 * __apple_build_version__, which upstream LLVM Clang does not.
 *
 * It must be detected separately because its __clang_major__ follows Apple's own numbering, which does not
 * correspond to upstream LLVM releases. Deciding whether a feature is available from __clang_major__ alone is
 * therefore wrong on this compiler.
 */
#if defined(__apple_build_version__)
#define CEU_COMPILER_IS_APPLE_CLANG
#endif

#if (defined(CEU_COMPILER_IS_APPLE_CLANG) && !defined(CEU_COMPILER_NAME))
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_APPLE_CLANG
#endif

#endif /* CEU_CC_APPLE_CLANG_H */
