#include "ceu_check/ceu_check_cc.h"
#include "ceu_check/ceu_check_c_std.h"
#include "ceu_check/ceu_check_c_stdlib.h"

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    {
        char* info = ceu_interpret_c_std_version();
        if (info == NULL) {
            fprintf(stderr, "Failed to get C standard information.\n");
            return EXIT_FAILURE;
        }
        printf("%s", info);
        free(info);
    }

    {
        char* info = ceu_check_get_compiler_info();
        if (info == NULL) {
            fprintf(stderr, "Failed to get compiler information.\n");
            return EXIT_FAILURE;
        }
        printf("%s", info);
        free(info);
    }
    {
        char* info = ceu_interpret_c_stdlib_version();
        if (info == NULL) {
            fprintf(stderr, "Failed to get C standard library information.\n");
            return EXIT_FAILURE;
        }
        printf("%s", info);
        free(info);
    }

    return EXIT_SUCCESS;
}
