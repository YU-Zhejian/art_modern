#include <stdlib.h> // abort()
#define TB_ASSERT(i) \
    if (!(i)) {      \
        abort();     \
    }

#define TEST_BOOL_F(BOOL_TYPE, INIT_TRUE, INIT_FALSE) \
    TB_ASSERT(sizeof(BOOL_TYPE) == 1)                 \
    BOOL_TYPE t = INIT_TRUE;                          \
    BOOL_TYPE f = INIT_FALSE;                         \
    TB_ASSERT(t);                                     \
    TB_ASSERT(t == 1);                                \
    TB_ASSERT(f == 0);                                \
    TB_ASSERT(!f);                                    \
    TB_ASSERT(!(!t));                                 \
    TB_ASSERT(!(t && f));                             \
    TB_ASSERT(t || f);                                \
    TB_ASSERT(t&& t);                                 \
    TB_ASSERT(!(f && f));                             \
    TB_ASSERT(!(f || f));\
