#include "ceu_check/ceu_check_c_stdlib.h"

#include "ceu_check/utils/c_snprintf.h"
#include "ceu_check/ceu_check_c_stdlib_macro.h"
#include "ceu_check/ceu_constants.h" /* NOLINT: for CEU_UNDEFINED */

#include <stdlib.h>

char* ceu_interpret_c_stdlib_version(void)
{
    char* buffer = calloc(1024, sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }

    char* ptr = buffer;
    int remaining = 1024;
    int written = 0;

    written = CEU_SNPRINTF(ptr, remaining, "Compile-time C stdlib: %s\n", CEU_LIBC_NAME);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#ifdef CEU_LIBC_IS_GNU
#ifdef __GLIBC__
    /* NOLINTNEXTLINE: __GLIBC_MINOR__ in <feature.h> should have been included if possible. */
    written = CEU_SNPRINTF(ptr, remaining, "\tglibc version number: %d.%d\n", __GLIBC__, __GLIBC_MINOR__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "\tglibc version number: %s\n", CEU_UNDEFINED);
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif
#ifdef CEU_LIBC_IS_UCLIBC
#if defined(__UCLIBC_MAJOR__) && defined(__UCLIBC_MINOR__) && defined(__UCLIBC_SUBLEVEL__)
    written = CEU_SNPRINTF(
        ptr, remaining, "\tuClibc version number: %d.%d.%d\n", __UCLIBC_MAJOR__, __UCLIBC_MINOR__, __UCLIBC_SUBLEVEL__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "\tuClibc version number: %s\n", CEU_UNDEFINED);
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif
#ifdef CEU_LIBC_IS_NEWLIB
#if defined(__NEWLIB__) && defined(__NEWLIB_MINOR__) && defined(__NEWLIB_PATCHLEVEL__)
    written = CEU_SNPRINTF(
        ptr, remaining, "\tNewlib version number: %d.%d.%d\n", __NEWLIB__, __NEWLIB_MINOR__, __NEWLIB_PATCHLEVEL__);
#else
    written = CEU_SNPRINTF(ptr, remaining, "\tNewlib version number: %s\n", CEU_UNDEFINED);
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }
#endif
    return buffer;
}
