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

int ceu_is_system_16_bit(void)
{
#if defined(CEU_ARCHITECTURE_16_BIT)
    return 1;
#elif defined(CEU_ARCHITECTURE_32_BIT) || defined(CEU_ARCHITECTURE_64_BIT)
    return 0;
#else
    /* Fallback method */
    return sizeof(void*) == 2;
#endif
}

int ceu_is_system_32_bit(void)
{
#if defined(CEU_ARCHITECTURE_32_BIT)
    return 1;
#elif defined(CEU_ARCHITECTURE_64_BIT) || defined(CEU_ARCHITECTURE_16_BIT)
    return 0;
#else
    /* Fallback method */
    return sizeof(void*) == 4;
#endif
}

int ceu_is_system_64_bit(void)
{
#if defined(CEU_ARCHITECTURE_64_BIT)
    return 1;
#elif defined(CEU_ARCHITECTURE_32_BIT) || defined(CEU_ARCHITECTURE_16_BIT)
    return 0;
#else
    /* Fallback method */
    return sizeof(void*) == 8;
#endif
}
