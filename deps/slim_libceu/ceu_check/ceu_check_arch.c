
#include "ceu_check/ceu_check_arch_macros.h"
#include "ceu_check/ceu_check_arch.h"


int ceu_is_big_endian(void)
{
#if ((defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)) || defined(__LITTLE_ENDIAN__)               \
    || defined(__LITTLE_ENDIAN) || defined(CEU_ARCHITECTURE_X86_64) || defined(CEU_ARCHITECTURE_I386)
    return 0;
#elif ((defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)) || defined(__BIG_ENDIAN__)                   \
    || defined(__BIG_ENDIAN)
    return 1;
#else /* Fallback method */
    unsigned int x = 1;
    return *((char*)&x) == 0;
#endif
}

int ceu_is_little_endian(void) /* NOLINT */ { return !ceu_is_big_endian(); /* NOLINT */ }

char* ceu_check_get_target_triplet_part_one(void) /* NOLINT */
{
#if defined(CEU_ARCHITECTURE_X86_64)
    return "x86_64";
#elif defined(CEU_ARCHITECTURE_I386)
    return "i386";
#endif
    return "unknown";
}
