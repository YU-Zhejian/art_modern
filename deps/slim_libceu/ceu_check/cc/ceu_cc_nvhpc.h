#ifndef CEU_CC_NVHPC_H
#define CEU_CC_NVHPC_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif

#if defined(__NVCOMPILER_MAJOR__)
#define CEU_COMPILER_IS_NVHPC
#define CEU_COMPILER_NAME "NVidia High Performance Computing (NVHPC) Toolkit NVC++"
#endif

#endif
