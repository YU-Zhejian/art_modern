#include "ceu_check/ceu_check_c_std.h"

#include "ceu_check/ceu_check_c_std_macro.h"
#include "ceu_check/utils/c_snprintf.h"

#include <stdlib.h>

char* ceu_interpret_c_std_version(void)
{
    char* buffer = calloc(1024, sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }

    char* ptr = buffer;
    int remaining = 1024;
    int written = 0;

    written = CEU_SNPRINTF(ptr, remaining, "C compile-time std.: ver. %s", CEU_C_STD);
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

#ifdef CEU_C_STD_VERSION_MACRO
    written = CEU_SNPRINTF(ptr, remaining, " (%ld)\n", (long)CEU_C_STD_VERSION_MACRO);
#else
    written = CEU_SNPRINTF(ptr, remaining, " (%s)\n", CEU_UNDEFINED);
#endif
    if (written > 0 && written < remaining) {
        ptr += written;
        remaining -= written;
    } else {
        free(buffer);
        return NULL;
    }

    return buffer;
}
