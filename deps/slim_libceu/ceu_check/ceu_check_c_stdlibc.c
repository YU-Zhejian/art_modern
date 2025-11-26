/**
 * TODO: Check libc
 */

#include <limits.h> /* NOLINT for __GLIBC__ */

#if defined(__GLIBC__) || defined(__GNU_LIBRARY__)
#define CEU_LIBC_IS_GNU
#else
#define CEU_LIBC_IS_UNKNOWN
#endif
