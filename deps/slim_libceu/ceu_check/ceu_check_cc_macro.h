/*!
 * @file ceu_check_cc_macro.h
 * @brief Header to get compiler information at compile time.
 * This file defines macros only.
 *
 * @see https://sourceforge.net/p/predef/wiki/Architectures/
 */

#ifndef CEU_CHECK_CC_MACRO_H
#define CEU_CHECK_CC_MACRO_H

#ifdef CEU_UNDER_DOXYGEN

/*!
 * @def CEU_COMPILER_IS_AOCC
 * @brief If this macro is defined, the compiler should support AMD Optimizing C/C++ and Fortran Compilers (AOCC)
 * features.
 * @see https://www.amd.com/en/developer/aocc.html
 * @see [here](https://www.techpowerup.com/233548/amd-ryzen-optimized-c-and-c-compilers-improve-performance) for a
 * benchmark comparing AOCC and other compilers.
 */
#define CEU_COMPILER_IS_AOCC

/*!
 * @def CEU_COMPILER_IS_LLVM
 * @brief If this macro is defined, the compiler should support LLVM features.
 */
#define CEU_COMPILER_IS_LLVM

/*!
 * @def CEU_COMPILER_IS_CLANG
 * @brief If this macro is defined, the compiler should support Clang features.
 */
#define CEU_COMPILER_IS_CLANG

/*!
 * @def CEU_COMPILER_IS_EDG
 * @brief If this macro is defined, the compiler should support Edison Design Group (EDG) features.
 * @see https://www.edg.com/
 */
#define CEU_COMPILER_IS_EDG

/*!
 * @def CEU_COMPILER_IS_GCC
 * @brief If this macro is defined, the compiler should support GCC/G++ features.
 */
#define CEU_COMPILER_IS_GCC

/*!
 *@def CEU_COMPILER_IS_ICC
 *@brief If this macro is defined, the compiler should support Intel(R) C++ Compiler Classic (ICC) features.
 */
#define CEU_COMPILER_IS_ICC

/*!
 *@def CEU_COMPILER_IS_INTEL_CLANG
 *@brief If this macro is defined, the compiler should support Intel(R) oneAPI DPC++/C++ Compiler (ICX) features.
 */
#define CEU_COMPILER_IS_INTEL_CLANG

/*!
 *@def CEU_COMPILER_IS_MSVC
 *@brief If this macro is defined, the compiler should support MSVC features.
 */
#define CEU_COMPILER_IS_MSVC

/*!
 * @def CEU_VISUAL_STUDIO_VER
 * @brief Human-readable Microsoft Visual Studio version.
 *
 * @see https://learn.microsoft.com/en-us/cpp/overview/compiler-versions?view=msvc-170
 * @see https://learn.microsoft.com/en-us/cpp/preprocessor/predefined-macros?view=msvc-170
 * @see https://devblogs.microsoft.com/cppblog/visual-c-compiler-version/
 */
#define CEU_VISUAL_STUDIO_VER /* implementation-defined */

/*!
 * @def CEU_COMPILER_IS_NVHPC
 * @brief If this macro is defined, the compiler should support NVidia HPC (NVHPC) Toolkit NVC++features.
 */
#define CEU_COMPILER_IS_NVHPC

/*!
 * @def CEU_COMPILER_IS_PGIC
 * @brief If this macro is defined, the compiler should support Portland Group (PGI) compiler features.
 */
#define CEU_COMPILER_IS_PGIC

/*!
 * @def CEU_COMPILER_IS_TCC
 * @brief If this macro is defined, the compiler should support Tiny C Compiler (TCC) features.
 */
#define CEU_COMPILER_IS_TCC

/*!
 * @def CEU_COMPILER_IS_BORLAND
 * @brief Indicator for Borland C/C++ compilers.
 */
#define CEU_COMPILER_IS_BORLAND

/*!
 * @def CEU_CPPB_VERSION
 * @brief Indicator for Borland/Embarcadero C++Builder version.
 *
 * @see https://docwiki.embarcadero.com/RADStudio/Alexandria/en/Compiler_Versions
 * @see https://docwiki.embarcadero.com/RADStudio/Alexandria/en/Predefined_Macros
 */
#define CEU_CPPB_VERSION /* implementation-defined */

/*!
 * @def CEU_COMPILER_NAME
 * @brief Human-readable compiler name.
 */
#define CEU_COMPILER_NAME /* implementation-defined */

#else

// Undefine conflicting macros
#ifdef CEU_COMPILER_NAME
#undef CEU_COMPILER_NAME
#endif

// Inclusion order of following lines should be preserved
// clang-format off
#include "ceu_check/cc/ceu_cc_turboc.h" // NOLINT
#include "ceu_check/cc/ceu_cc_icc.h" // NOLINT
#include "ceu_check/cc/ceu_cc_intel_clang.h" // NOLINT
#include "ceu_check/cc/ceu_cc_amd_clang.h" // NOLINT
#include "ceu_check/cc/ceu_cc_msvc.h" // NOLINT
#include "ceu_check/cc/ceu_cc_nvhpc.h" // NOLINT
#include "ceu_check/cc/ceu_cc_pgic.h" // NOLINT
#include "ceu_check/cc/ceu_cc_tcc.h" // NOLINT
#include "ceu_check/cc/ceu_cc_clang.h" // NOLINT
#include "ceu_check/cc/ceu_cc_gcc.h" // NOLINT
#include "ceu_check/cc/ceu_cc_edg.h" // NOLINT
#include "ceu_check/cc/ceu_cc_unknown.h" // NOLINT
// clang-format on
#endif
#endif // CEU_CHECK_CC_MACRO_H
