#include "ceu_check/ceu_check_c_std.h"
#include "ceu_check/c_snprintf.h"
#include "ceu_check/ceu_check_c_cxx_std_macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* ceu_interpret_c_std_version()
{
    char* buffer = calloc(1024, sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }

    char cstd_macro[256];
    char* ptr = buffer;
    int remaining = 1024;
    int written;

    written = CEU_SNPRINTF(ptr, remaining, "Compile-time C std.: ver. %s", CEU_C_STD);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

#ifdef CEU_C_STD_VERSION_MACRO
    CEU_SNPRINTF(cstd_macro, sizeof(cstd_macro), "%ld", (long)CEU_C_STD_VERSION_MACRO);
#else
    CEU_SNPRINTF(cstd_macro, sizeof(cstd_macro), "__STDC_VERSION__ undefined");
#endif

    written = CEU_SNPRINTF(ptr, remaining, " (%s)\n", cstd_macro);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

    return buffer;
}
