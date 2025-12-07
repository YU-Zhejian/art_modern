/**
 *
 * TODO: Support more architectures.
 */
#ifndef CEU_CHECK_ARCH_H
#define CEU_CHECK_ARCH_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Detect the endianness of the target machine.
 *
 * @return 1 if the machine is big endian, 0 otherwise.
 */
int ceu_is_big_endian(void);

/**
 * @return Reverse of ceu_is_big_endian().
 */
int ceu_is_little_endian(void);

/**
 * Get the target architecture triplet's first part.
 *
 * @return A string representing the architecture, or "unknown" if not detected.
 */
char* ceu_check_get_target_triplet_part_one(void) ;

#ifdef __cplusplus
}
#endif
#endif