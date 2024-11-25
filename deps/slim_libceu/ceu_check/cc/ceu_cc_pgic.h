// Include this file after <ceu_cc_nvhpc.h>!
#ifndef CEU_CC_PGIC_H
#define CEU_CC_PGIC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif

#if defined(__PGI)
#define CEU_COMPILER_IS_PGIC
#if !defined(CEU_COMPILER_NAME)
#define CEU_COMPILER_NAME "Portland Group (PGI) C/C++ compilers"
#endif
#endif

#endif
