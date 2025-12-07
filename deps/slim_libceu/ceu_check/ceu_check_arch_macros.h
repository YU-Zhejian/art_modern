#ifndef CEU_CHECK_ARCH_MACRO_H
#define CEU_CHECK_ARCH_MACRO_H

#include "ceu_check/ceu_check_c_cxx_compiler_macro.h" /* NOLINT: for CEU_COMPILER_IS_TI */
#include "libceu_stddef.h" /* NOLINT: CEU_HAVE_INCLUDE_ENDIAN_H */

#if defined(CEU_HAVE_INCLUDE_ENDIAN_H) && (CEU_HAVE_INCLUDE_ENDIAN_H == 1)
#include <endian.h>
#if (!defined(BYTE_ORDER)) || (!defined(LITTLE_ENDIAN)) || (!defined(BIG_ENDIAN)) || (!defined(PDP_ENDIAN))
#error "<endian.h> is included but BYTE_ORDER, LITTLE_ENDIAN, BIG_ENDIAN or PDP_ENDIAN is not defined."
#endif
#define CEU_INCLUDED_ENDIAN_H
#endif

/* GCC 4.6+ and Clang provide __BYTE_ORDER__ */
#ifdef CEU_COMPILER_IS_GCC
#if (defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__) && defined(__ORDER_PDP_ENDIAN__))
#define CEU_GCC_LIKE_ENDIAN_MACROS
#endif
#endif

/* ---------------------------------------------------------------------- */
/** Define architecture macros
 */

#if defined(__x86_64__) || defined(__x86_64) || defined(_M_AMD64) || defined(__amd64__) || defined(__amd64)
#define CEU_ARCHITECTURE_X86_64
#elif defined(i386) || defined(__i386) || defined(__i386__) || defined(__i486__) || defined(__i586__)                  \
    || defined(__i686__) || defined(_M_IX86) || defined(__X86__) || defined(_X86_)
#define CEU_ARCHITECTURE_I386
#else
#define CEU_ARCHITECTURE_UNKNOWN
#endif

/* ---------------------------------------------------------------------- */
/** Define endianness macros
 */

/* Firstly, we would error on PDP-endian, which is not supported here */
#if (defined(CEU_INCLUDED_ENDIAN_H) && (BYTE_ORDER == PDP_ENDIAN))                                                     \
    || (defined(CEU_GCC_LIKE_ENDIAN_MACROS) && (__BYTE_ORDER__ == __ORDER_PDP_ENDIAN__))
#error "PDP-endian is not supported."
#endif

/**
 * Define CEU_COMPILE_TIME_IS_BIG_ENDIAN if the target machine is detected to be big-endian at compile-time.
 */
#ifdef CEU_COMPILE_TIME_IS_BIG_ENDIAN
#undef CEU_COMPILE_TIME_IS_BIG_ENDIAN
#endif
#ifdef CEU_COMPILE_TIME_IS_LITTLE_ENDIAN
#undef CEU_COMPILE_TIME_IS_LITTLE_ENDIAN
#endif
/* clang-format off */
/**
 * Note that on some machines. __FLOAT_WORD_ORDER__ and __BYTE_ORDER__ may differ.
 *
 * ARM Embedded: see <https://developer.arm.com/documentation/dui0472/m/Compiler-specific-Features/Predefined-macros?lang=en>
 */
#if \
    0 /* Placeholder for future LE architectures */                                                                    \
    || /* Known LE Arch */ defined(CEU_ARCHITECTURE_X86_64)                                                            \
    || /* Known LE Arch */ defined(CEU_ARCHITECTURE_I386)                                                              \
    || /* <endian.h> */ (defined(CEU_INCLUDED_ENDIAN_H) && (BYTE_ORDER == LITTLE_ENDIAN))                              \
    || /* GCC */  (defined(CEU_GCC_LIKE_ENDIAN_MACROS) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__))                 \
    || /* TI ARMCL */ (defined(CEU_COMPILER_IS_TI) && defined(__LITTLE_ENDIAN__) && (__LITTLE_ENDIAN__ == 1))          \
    || /* Clang */ (defined(CEU_COMPILER_IS_CLANG) && defined(__LITTLE_ENDIAN__) && (__LITTLE_ENDIAN__ == 1))          \
    || /* TI ARMCL */ (defined(CEU_COMPILER_IS_TI) && defined(__little_endian__) && (__little_endian__ == 1))          \
    || /* TI CL6X */ (defined(CEU_COMPILER_IS_TI) && defined(_LITTLE_ENDIAN) && (_LITTLE_ENDIAN == 1))                 \
    || /* TIARMCLANG */ (defined(CEU_COMPILER_IS_TI) && (!defined(__ARM_BIG_ENDIAN)))                                  \
    || /* ARM Embedded */ (defined(CEU_COMPILER_IS_ARM_COMPILER_EMBEDDED) && (!defined(__BIG_ENDIAN)))                 \
    || 0
#define CEU_COMPILE_TIME_IS_LITTLE_ENDIAN
#elif \
    0 /* Placeholder for future BE architectures */                                                                    \
    || /* <endian.h> */ (defined(CEU_INCLUDED_ENDIAN_H) && (BYTE_ORDER == BIG_ENDIAN))                                 \
    || /* GCC */  (defined(CEU_GCC_LIKE_ENDIAN_MACROS) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__))                    \
    || /* TI ARMCL */ (defined(CEU_COMPILER_IS_TI) && defined(__BIG_ENDIAN__) && (__BIG_ENDIAN__ == 1))                \
    || /* Clang */ (defined(CEU_COMPILER_IS_CLANG) && defined(__BIG_ENDIAN__) && (__BIG_ENDIAN__ == 1))                \
    || /* TI ARMCL */ (defined(CEU_COMPILER_IS_TI) && defined(__big_endian__) && (__big_endian__ == 1))                \
    || /* TIARMCLANG */ (defined(CEU_COMPILER_IS_TI) && defined(__ARM_BIG_ENDIAN))                                     \
    || /* ARM Embedded */ (defined(CEU_COMPILER_IS_ARM_COMPILER_EMBEDDED) && defined(__BIG_ENDIAN))                    \
    || 0
#define CEU_COMPILE_TIME_IS_BIG_ENDIAN
#endif
/* clang-format on */

#if (!defined(CEU_COMPILE_TIME_IS_LITTLE_ENDIAN)) && (!defined(CEU_COMPILE_TIME_IS_BIG_ENDIAN))
#define CEU_COMPILE_TIME_IS_UNKNOWN_ENDIAN
#endif

#ifdef CEU_INCLUDED_ENDIAN_H
#undef CEU_INCLUDED_ENDIAN_H
#endif

#endif /* CEU_CHECK_ARCH_MACRO_H */
