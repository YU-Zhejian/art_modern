
#if defined(__x86_64__) || defined(__x86_64) || defined(_M_AMD64) || defined(__amd64__) || defined(__amd64)
#define CEU_ARCHITECTURE_X86_64
#elif defined(i386) || defined(__i386) || defined(__i386__) || defined(__i486__) || defined(__i586__)                  \
    || defined(__i686__) || defined(_M_IX86) || defined(__X86__) || defined(_X86_)
#define CEU_ARCHITECTURE_I386
#endif

/**
 * Detect endianness of the target machine.
 *
 * @return 1 if the machine is big endian, 0 otherwise.
 */
int ceu_is_big_endian(void)
{
#if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) || defined(__LITTLE_ENDIAN__)               \
    || defined(__LITTLE_ENDIAN) || defined(CEU_ARCHITECTURE_X86_64) || defined(CEU_ARCHITECTURE_I386)
    return 0;
#elif (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__) || defined(__BIG_ENDIAN__)                   \
    || defined(__BIG_ENDIAN)
    return 1;
#else /* Fallback method */
    unsigned int x = 1;
    return *((char*)&x) == 0;
#endif
}

/**
 * @return Reverse of ceu_is_big_endian().
 */
int ceu_is_little_endian(void) { return !ceu_is_big_endian(); }

/**
 * Get the target architecture triplet's first part.
 *
 * @return A string representing the architecture, or "unknown" if not detected.
 */
char* ceu_check_get_target_triplet_part_one(void)
{
#if defined(CEU_ARCHITECTURE_X86_64)
    return "x86_64";
#elif defined(CEU_ARCHITECTURE_I386)
    return "i386";
#endif
    return "unknown";
}
