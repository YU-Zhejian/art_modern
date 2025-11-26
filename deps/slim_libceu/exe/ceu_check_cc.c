#include "ceu_check/ceu_check_cc.h"
#include "libceu_stddef.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    char* info = ceu_check_get_compiler_info();
    if (info != NULL) {
        printf("%s", info);
        free(info);
    } else {
        fprintf(stderr, "Failed to get compiler information.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}