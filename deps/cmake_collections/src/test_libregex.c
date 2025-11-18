#include <regex.h>
#include <stddef.h>
#include <stdlib.h>

const char* test_string = "this is a test string";
const char* pattern = "test.*";
const char* non_matching_string = "no match here";

int main(void)
{
    regex_t re;
    int ret = regcomp(&re, "test.*", REG_EXTENDED);
    if (ret != 0) {
        return EXIT_FAILURE;
    }
    // Test the regex
    if (regexec(&re, test_string, 0, NULL, 0) != 0) {
        regfree(&re);
        return EXIT_FAILURE; // Should match
    }
    if (regexec(&re, non_matching_string, 0, NULL, 0) == 0) {
        regfree(&re);
        return EXIT_FAILURE; // Should not match
    }

    regfree(&re);
    return EXIT_SUCCESS;
}
