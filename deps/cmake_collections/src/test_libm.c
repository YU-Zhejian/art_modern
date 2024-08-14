#include <math.h>

int main(void)
{
    volatile int s = (int)sqrt(1);
    return s - 1;
}
