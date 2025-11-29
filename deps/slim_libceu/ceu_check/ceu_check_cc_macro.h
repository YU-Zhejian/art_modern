/*!
 * @file ceu_check_cc_macro.h
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

#include "ceu_check/ceu_compiler_names.h"

/* The inclusion order of the following lines should be preserved to avoid conflicts. */
/* clang-format off */
#include "ceu_check/cc/ceu_cc_turboc.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_icc.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_arm_compiler_linux.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_arm_compiler_embedded.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_intel_clang.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_amd_clang.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_msvc.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_nvhpc.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_pgic.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_tcc.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_clang.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_gcc.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_edg.h" /* NOLINT */
#include "ceu_check/cc/ceu_cc_unknown.h" /* NOLINT */
/* clang-format on */

#endif /* CEU_CHECK_CC_MACRO_H */
