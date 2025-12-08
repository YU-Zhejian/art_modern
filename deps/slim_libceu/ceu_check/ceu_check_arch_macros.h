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
#if (defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__)                      \
    && defined(__ORDER_PDP_ENDIAN__))
#define CEU_GCC_LIKE_ENDIAN_MACROS
#endif
#endif




/* ---------------------------------------------------------------------- */
/** Define architecture macros
 */

#ifdef CEU_ARCH_NAME
#undef CEU_ARCH_NAME
#endif


/* clang-format off */

#if defined(__x86_64__) || defined(__x86_64) || defined(_M_AMD64) || defined(__amd64__) || defined(__amd64)
#define CEU_ARCHITECTURE_X86_64
#define CEU_ARCH_NAME "x86_64"
#elif 0                                                                                                                \
    || /* i386 */ defined(i386) || defined(__i386) || defined(__i386__)                                                \
    || /* i?86 */ defined(__i486__) || defined(__i586__) || defined(__i686) || defined(__i686__)                       \
    || /* Other definitions */ defined(_M_IX86) || defined(__X86__) || defined(_X86_)
#define CEU_ARCHITECTURE_I386
#define CEU_ARCH_NAME "i386"
#elif defined(__aarch64__) || defined(_M_ARM64)
#define CEU_ARCHITECTURE_AARCH64
#define CEU_ARCH_NAME "aarch64"
#elif 0                                                                                                                \
    || defined(__arm) || defined(__arm__) || defined(_M_ARM) || defined(__TI_ARM__)                                    \
    || /* ARM arch <= 7 is 32-bit */ (defined(__ARM_ARCH) && __ARM_ARCH <= 7)                                          \
    || /* ARM arch <= 7 is 32-bit */ (defined(__ARM_ARCH__) && __ARM_ARCH__ <= 7)                                      \
    || /* ARM v7 */ defined(__ARM_ARCH_7__) || defined(__ARM_ARCH_7A__) || defined(__ARM_ARCH_7R__)                    \
    || /* ARM v7 */ defined(__ARM_ARCH_7M__) || defined(__TI_ARM_V7R5__) || defined(__TI_ARM_V7R4__)                   \
    || /* ARM v7 */ defined(__TI_ARM_V7M4__) || defined(__TI_ARM_V7M3__) || defined(__TI_ARM_V7__)                     \
    || /* ARM v7 */ defined(__TI_ARM_V7A8__) ||  defined(__ARM_ARCH_7S__)                                              \
    || /* ARM v6 */ defined(__ARM_ARCH_6__) || defined(__ARM_ARCH_6J__) || defined(__ARM_ARCH_6ZK__)                   \
    || /* ARM v6 */ defined(__ARM_ARCH_6K__) || defined(__ARM_ARCH_6Z__) || defined(__ARM_ARCH_6KZ__)                  \
    || /* ARM v6 */ defined(__ARM_ARCH_6T2__) || defined(__TI_ARM_V6M0__) || defined(__TI_ARM_V6__)                    \
    || /* ARM v5 */ defined(__ARM_ARCH_5TE__) || defined(__ARM_ARCH_5TEJ__) || defined(__TI_ARM_V5__)                  \
    || /* ARM v5 */ defined(__ARM_ARCH_5__) || defined(__ARM_ARCH_5E__) || defined(__ARM_ARCH_5T__)                    \
    || /* ARM v4 */ defined(__ARM_ARCH_4T__) || defined(__ARM_ARCH_4__) ||  || defined(__TI_ARM_V4__)                  \
    || /* ARM v4 */ defined(__TARGET_ARM_4T)                                                                           \
    || /* ARM v3 */ defined(__ARM_ARCH_3__) || defined(__ARM_ARCH_3M__)                                                \
    || /* ARM v2 */ defined(__ARM_ARCH_2__)                                                                            \
    || 0
#define CEU_ARCHITECTURE_ARM
#define CEU_ARCH_NAME "arm"
#elif (defined(__loongarch_grlen) && (__loongarch_grlen == 64)) || defined(__loongarch64)
#define CEU_ARCHITECTURE_LOONGARCH64
#define CEU_ARCH_NAME "loongarch64"
#elif defined(__loongarch_grlen) && (__loongarch_grlen == 32)
#define CEU_ARCHITECTURE_LOONGARCH32
#define CEU_ARCH_NAME "loongarch32"
#elif defined(__riscv_xlen) && (__riscv_xlen == 64)
#define CEU_ARCHITECTURE_RISCV64
#define CEU_ARCH_NAME "riscv64"
#elif defined(__riscv_xlen) && (__riscv_xlen == 32)
#define CEU_ARCHITECTURE_RISCV32
#define CEU_ARCH_NAME "riscv32"
#else
/** TODO: systemz, sparc, sparcel, sparcv9, ppc32, ppc32le, ppc64, ppc64le, mips, mipsel, mips64, mips64el
 * m68k */
#define CEU_ARCHITECTURE_UNKNOWN
#define CEU_ARCH_NAME "unknown"
#endif


#if 0                                                                                                                  \
    || /* Known 64-bit */ defined(__ARM_64BIT_STATE)                                                                   \
    || /* Known 64-bit */ defined(CEU_ARCHITECTURE_AARCH64)                                                            \
    || /* Known 64-bit */ defined(CEU_ARCHITECTURE_X86_64)                                                             \
    || /* Known 64-bit */ defined(CEU_ARCHITECTURE_LOONGARCH64)                                                        \
    || /* Known 64-bit */ defined(CEU_ARCHITECTURE_RISCV64)                                                            \
    || 0
#define CEU_ARCHITECTURE_64_BIT
#elif 0                                                                                                                \
    || /* Known 32-bit */ defined(__ARM_32BIT_STATE)                                                                   \
    || /* Known 32-bit */ defined(CEU_ARCHITECTURE_ARM)                                                                \
    || /* Known 32-bit */ defined(CEU_ARCHITECTURE_I386)                                                               \
    || /* Known 32-bit */ defined(CEU_ARCHITECTURE_LOONGARCH32)                                                        \
    || /* Known 32-bit */ defined(CEU_ARCHITECTURE_RISCV32)                                                            \
    || 0
#define CEU_ARCHITECTURE_32_BIT
#endif
/* clang-format on */

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
 *
 * ARM, RISC-V and Loongarch can be either little-endian or big-endian.
 */
#if \
    0 /* Placeholder for future LE architectures */                                                                    \
    || /* Known LE Arch */ defined(CEU_ARCHITECTURE_X86_64)                                                            \
    || /* Known LE Arch */ defined(CEU_ARCHITECTURE_I386)                                                              \
    || /* Known LE Arch */ defined(__ARMEL__)                                                                          \
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
    || /* Known BE Arch */ defined(__ARMEB__)                                                                          \
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

/* Clean up temporary macros */
#ifdef CEU_INCLUDED_ENDIAN_H
#undef CEU_INCLUDED_ENDIAN_H
#endif
#ifdef CEU_GCC_LIKE_ENDIAN_MACROS
#undef CEU_GCC_LIKE_ENDIAN_MACROS
#endif

#endif /* CEU_CHECK_ARCH_MACRO_H */
