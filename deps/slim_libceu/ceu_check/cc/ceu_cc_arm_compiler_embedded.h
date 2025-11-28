/**
 * Include before EDG
 *
 * NOTE: This compiler is propertairy software from Arm Ltd. We did not purchase a copy and the following
 * information is based on publicly available documentation.
 *
 */

#ifndef CEU_CC_ARM_COMPILER_EMBEDDED_H
#define CEU_CC_ARM_COMPILER_EMBEDDED_H

#ifndef CEU_CHECK_CC_MACRO_H
#error "Do not include this file, include <ceu_check/ceu_check_cc_macro.h> instead!"
#endif
#if defined(__ARMCC_VERSION) || defined(__CC_ARM)
#define CEU_COMPILER_IS_ARM_COMPILER_EMBEDDED
#endif
#if (defined(CEU_COMPILER_IS_ARM_COMPILER_EMBEDDED) && !defined(CEU_COMPILER_NAME))
#define CEU_COMPILER_NAME "Arm Compiler for Embedded (ACOMPE)"
#endif
#endif /* CEU_CC_ARM_COMPILER_EMBEDDED_H */
