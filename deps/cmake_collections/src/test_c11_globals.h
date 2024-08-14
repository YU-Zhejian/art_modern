#include "test_c99_globals.h"
/* Does the compiler advertise C11 conformance? */
#if !defined __STDC_VERSION__ || __STDC_VERSION__ < 201112L
#error "Compiler does not advertise C11 conformance"
#endif

// Check _Alignas.
char _Alignas(double) aligned_as_double;
char _Alignas(0) no_special_alignment;
extern char aligned_as_int;
char _Alignas(0) _Alignas(int) aligned_as_int;

// Check _Alignof.
enum {
    int_alignment = _Alignof(int),
    int_array_alignment = _Alignof(int[100]),
    char_alignment = _Alignof(char)
};
_Static_assert(0 < -_Alignof(int), "_Alignof is signed");

// Check _Noreturn.
int _Noreturn does_not_return(void)
{
    for (;;)
        continue;
}

// Check _Static_assert.
struct test_static_assert {
    int x;
    _Static_assert(sizeof(int) <= sizeof(long int),
        "_Static_assert does not work in struct");
    long int y;
};

// Check UTF-8 literals.
#define u8 syntax error !
char const utf8_literal[] = u8"happens to be ASCII"
                            "another string";

// Check duplicate typedefs.
typedef long* long_ptr;
typedef long int* long_ptr;
typedef long_ptr long_ptr;

// Anonymous structures and unions -- taken from C11 6.7.2.1 Example 1.
struct anonymous {
    union {
        struct {
            int i;
            int j;
        };
        struct {
            int k;
            long int l;
        } w;
    };
    int m;
} v1;

int c11main(void)
{
    int ok = 0;
    _Static_assert((offsetof(struct anonymous, i)
                       == offsetof(struct anonymous, w.k)),
        "Anonymous union alignment botch");
    v1.i = 2;
    v1.w.k = 5;
    ok |= v1.i != 5;
    return ok;
}
