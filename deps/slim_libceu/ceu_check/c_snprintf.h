#ifndef CEU_CHECK_C_SNPRINTF_H
#define CEU_CHECK_C_SNPRINTF_H
#ifdef CEU_COMPILER_IS_MSVC

/**  Enable MSVC secure CRT functions */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#endif
#include <stdio.h>

#ifdef CEU_COMPILER_IS_MSVC
/** Macro for _snprintf_s with truncation behavior similar to snprintf  */
#define CEU_SNPRINTF(buf, buf_size, ...) \
    _snprintf_s(buf, buf_size, _TRUNCATE, __VA_ARGS__)
#else
/** Macro for snprintf */
#define CEU_SNPRINTF(buf, buf_size, ...) \
    snprintf(buf, buf_size, __VA_ARGS__)
#endif

#endif /* CEU_CHECK_C_SNPRINTF_H */
