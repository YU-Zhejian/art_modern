#include <regex.h>

int main(void)
{
    regex_t re;
    int ret = regcomp(&re, "test", REG_EXTENDED);
    regfree(&re);
    return ret;
}