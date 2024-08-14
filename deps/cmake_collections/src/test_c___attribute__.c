/*!
 * @file test_c___attribute__.c
 * @author YU Zhejian
 * @brief Test whether __attribute__ is supported. Will test __attribute__((noreturn)) and __attribute__((unused)).
 * @version 0.1
 * @date 2024-04-28
 */
#include <stdlib.h>

__attribute__((noreturn)) static int f(int i)
{
    exit(0);
}

__attribute__((unused)) static int g(int i)
{
    return 1;
}

int main(void)
{
    f(1);
    return 0;
}
