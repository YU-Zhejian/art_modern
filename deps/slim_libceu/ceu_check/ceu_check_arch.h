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

int ceu_is_system_16_bit(void);
int ceu_is_system_32_bit(void);
int ceu_is_system_64_bit(void);

#ifdef __cplusplus
}
#endif
#endif
