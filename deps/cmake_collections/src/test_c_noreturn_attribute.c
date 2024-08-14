#include <stdlib.h>

[[noreturn]] static int f(int i)
{
    exit(0);
}

int main(void)
{
    f(1);
    return 0;
}
