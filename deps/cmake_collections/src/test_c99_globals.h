#include "test_c89_globals.h"
/* Does the compiler advertise C99 conformance? */
#if !defined __STDC_VERSION__ || __STDC_VERSION__ < 199901L
#error "Compiler does not advertise C99 conformance"
#endif

// See if C++-style comments work.

#include <stdbool.h>
extern int puts(const char*);
extern int printf(const char*, ...);
// extern int dprintf(int, const char*, ...);
// FIXME: Not working under Windows
extern void* malloc(size_t);
extern void free(void*);

// Check varargs macros.  These examples are taken from C99 6.10.3.5.
// dprintf is used instead of fprintf to avoid needing to declare
// FILE and stderr.
// #define debug(...) dprintf(2, __VA_ARGS__)
#define debug(...) // FIXME: See above
#define showlist(...) puts(#__VA_ARGS__)
#define report(test, ...) ((test) ? puts(#test) : printf(__VA_ARGS__))
static void
test_varargs_macros(void)
{
    int x = 1234;
    int y = 5678;
    debug("Flag");
    debug("X = %d\n", x);
    showlist(The first, second, and third items.);
    report(x > y, "x is %d but y is %d", x, y);
}

// Check long long types.
#define BIG64 18446744073709551615ull
#define BIG32 4294967295ul
#define BIG_OK (BIG64 / BIG32 == 4294967297ull && BIG64 % BIG32 == 0)
#if !BIG_OK
#error "your preprocessor is broken"
#endif
#if BIG_OK
#else
#error "your preprocessor is broken"
#endif
static long long int bignum = -9223372036854775807LL;
static unsigned long long int ubignum = BIG64;

struct incomplete_array {
    int datasize;
    double data[];
};

struct named_init {
    int number;
    const wchar_t* name;
    double average;
};

typedef const char* ccp;

static inline int
test_restrict(ccp restrict text)
{
    // Iterate through items via the restricted pointer.
    // Also check for declarations in for loops.
    for (unsigned int i = 0; *(text + i) != '\0'; ++i)
        continue;
    return 0;
}

// Check varargs and va_copy.
static bool
test_varargs(const char* format, ...)
{
    va_list args;
    va_start(args, format);
    va_list args_copy;
    va_copy(args_copy, args);

    const char* str = "";
    int number = 0;
    float fnumber = 0;

    while (*format) {
        switch (*format++) {
        case 's': // string
            str = va_arg(args_copy, const char*);
            break;
        case 'd': // int
            number = va_arg(args_copy, int);
            break;
        case 'f': // float
            fnumber = va_arg(args_copy, double);
            break;
        default:
            break;
        }
    }
    va_end(args_copy);
    va_end(args);

    return *str && number && fnumber;
}

int c99main(int argc, char** argv)
{
    int ok = 0;
    // Check bool.
    _Bool success = false;
    success |= (argc != 0);

    // Check restrict.
    if (test_restrict("String literal") == 0)
        success = true;
    char* restrict newvar = "Another string";

    // Check varargs.
    success &= test_varargs("s, d'\'' f .", "string", 65, 34.234);
    test_varargs_macros();

    // Check flexible array members.
    struct incomplete_array* ia = malloc(sizeof(struct incomplete_array) + (sizeof(double) * 10));
    ia->datasize = 10;
    for (int i = 0; i < ia->datasize; ++i)
        ia->data[i] = i * 1.234;
    // Work around memory leak warnings.
    free(ia);

    // Check named initializers.
    struct named_init ni = {
        .number = 34,
        .name = L"Test wide string",
        .average = 543.34343,
    };

    ni.number = 58;

    int dynamic_array[ni.number];
    dynamic_array[0] = argv[0][0];
    dynamic_array[ni.number - 1] = 543;

    // work around unused variable warnings
    ok |= (!success || bignum == 0LL || ubignum == 0uLL || newvar[0] == 'x' || dynamic_array[ni.number - 1] != 543);
    return ok;
}
