#ifndef CEU_INCLUDES_X86_SIMD_H
#define CEU_INCLUDES_X86_SIMD_H

#include "libceu_stddef.h" /* NOLINT: for CEU_HAVE_INCLUDE_INTRIN_H and CEU_HAVE_INCLUDE_X86INTRIN_H */

#include "ceu_check/ceu_check_arch_macros.h" /* NOLINT: for CEU_ARCHITECTURE_X86_64 and CEU_ARCHITECTURE_I386 */

#if defined(CEU_ARCHITECTURE_X86_64) || defined(CEU_ARCHITECTURE_I386)
#if defined(CEU_COMPILER_IS_MSVC)
#if defined(CEU_HAVE_INCLUDE_INTRIN_H) && (CEU_HAVE_INCLUDE_INTRIN_H == 1)
#include <intrin.h> /* NOLINT */
#endif
#elif defined(CEU_COMPILER_IS_GCC)
#if defined(CEU_HAVE_INCLUDE_X86INTRIN_H) && (CEU_HAVE_INCLUDE_X86INTRIN_H == 1)
#include <x86intrin.h> /* NOLINT */
#endif
#else
/* Unknown compiler */
#endif
#else
/* No SIMD support for non-x86_64 architectures yet */
#endif

#endif /* CEU_INCLUDES_X86_SIMD_H */
