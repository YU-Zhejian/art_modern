#ifndef TEST_AM_ULTRA_DEEP_D_TYPES_H
#define TEST_AM_ULTRA_DEEP_D_TYPES_H

#include <stdint.h>
#include <stdlib.h>

/** All hash must be 64-bit **/
typedef uint64_t hash_type;
/** All hash functions must return a 64-bit hash */
typedef hash_type (*hash_function)(const void* data, const size_t len);
#endif // TEST_AM_ULTRA_DEEP_D_TYPES_H
