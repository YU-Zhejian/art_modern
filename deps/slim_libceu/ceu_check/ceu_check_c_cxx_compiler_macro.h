/*!
 * @file ceu_check_c_cxx_compiler_macro.h
 * @brief Header to get compiler information at compile time.
 * This file defines macros only.
 *
 * @see https://sourceforge.net/p/predef/wiki/Architectures/
 */

#ifndef CEU_CHECK_CC_MACRO_H
#define CEU_CHECK_CC_MACRO_H

/* Undefine conflicting macros first */
#ifdef CEU_COMPILER_NAME
#undef CEU_COMPILER_NAME
#endif

#include "ceu_utils/ceu_compiler_names.h" /* NOLINT: for diverse compiler names used in the following includes */

/* The inclusion order of the following lines should be preserved to avoid conflicts. */
/* clang-format off */
#include "ceu_check/cc/01_turboc.h" /* NOLINT */
#include "ceu_check/cc/02_icc.h" /* NOLINT */
#include "ceu_check/cc/03_arm_compiler_linux.h" /* NOLINT */
#include "ceu_check/cc/04_arm_compiler_embedded.h" /* NOLINT */
#include "ceu_check/cc/05_intel_clang.h" /* NOLINT */
#include "ceu_check/cc/06_amd_clang.h" /* NOLINT */
#include "ceu_check/cc/07_msvc.h" /* NOLINT */
#include "ceu_check/cc/08_nvhpc.h" /* NOLINT */
#include "ceu_check/cc/09_pgic.h" /* NOLINT */
#include "ceu_check/cc/10_tcc.h" /* NOLINT */
#include "ceu_check/cc/11_ti.h" /* NOLINT */
#include "ceu_check/cc/30_clang.h" /* NOLINT */
#include "ceu_check/cc/40_gcc.h" /* NOLINT */
#include "ceu_check/cc/50_edg.h" /* NOLINT */
#include "ceu_check/cc/99_unknown.h" /* NOLINT */
/* clang-format on */

#endif /* CEU_CHECK_CC_MACRO_H */
