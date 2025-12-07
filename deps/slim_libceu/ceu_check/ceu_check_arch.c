#include "ceu_check/ceu_check_arch.h"
#include "ceu_check/ceu_check_arch_macros.h"

int ceu_is_big_endian(void)
{

#if defined(CEU_COMPILE_TIME_IS_LITTLE_ENDIAN)
    return 0;
#elif defined(CEU_COMPILE_TIME_IS_BIG_ENDIAN)
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
