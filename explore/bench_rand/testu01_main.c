#include <TestU01.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
    const char* temp_file = "/dev/stdin";
    char* filename = (char*)malloc(strlen(temp_file) + 1);
    if (filename == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }
    strncpy(filename, temp_file, strlen(temp_file) + 1);

    unif01_Gen* gen = ufile_CreateReadBin(filename, 4096);
    free(filename);

    swrite_Basic = FALSE;
    swrite_Host = FALSE;
    swrite_Parameters = FALSE;
    swrite_Collectors = FALSE;
    swrite_Classes = FALSE;
    swrite_Counters = FALSE;

    bbattery_SmallCrush(gen);

    swrite_Basic = FALSE;
    swrite_Host = FALSE;
    swrite_Parameters = FALSE;
    swrite_Collectors = FALSE;
    swrite_Classes = FALSE;
    swrite_Counters = FALSE;
    bbattery_pseudoDIEHARD(gen);

    unif01_DeleteExternGenBits(gen);

    return EXIT_SUCCESS;
}
