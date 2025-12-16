#ifndef CEU_CC_TCC_H
#define CEU_CC_TCC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_c_cxx_compiler_macro.h> instead!"
#endif

#if defined(__TINYC__)
#define CEU_COMPILER_IS_TCC
#define CEU_COMPILER_NAME CEU_COMPILER_NAME_TCC
#endif
#endif /* CEU_CC_TCC_H */
