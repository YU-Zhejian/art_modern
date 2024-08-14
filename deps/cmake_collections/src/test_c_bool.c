#include "test_bool.h"

#if defined(C_UNDERLINE_BOOL_KEYWORD)
#define CEU_TB_BOOL _Bool
#define CEU_TB_T ((_Bool) 1)
#define CEU_TB_F ((_Bool) 0)

#elif defined(C_STDBOOL)
#include <stdbool.h>
#ifndef __bool_true_false_are_defined
#error "__bool_true_false_are_defined not defined!"
#endif

#if __bool_true_false_are_defined != 1
#error "__bool_true_false_are_defined != 1!"
#endif

#define CEU_TB_BOOL bool
#define CEU_TB_T true
#define CEU_TB_F false

#else
#define CEU_TB_BOOL int
#define CEU_TB_T 1
#define CEU_TB_F 0

#endif

void test1()
{
    TEST_BOOL_F(CEU_TB_BOOL, CEU_TB_T, CEU_TB_F);
}
void test2()
{
    TEST_BOOL_F(CEU_TB_BOOL, 1, 0);
}
void test3()
{
    TEST_BOOL_F(CEU_TB_BOOL, 4, 0);
}
void test4()
{
    TEST_BOOL_F(CEU_TB_BOOL, -4, 0);
}

void test5()
{
    TEST_BOOL_F(CEU_TB_BOOL, 0.5, 0);
}

int main(void)
{
    test2();
    test3();
    test4();
    test5();
    return 0;
}
