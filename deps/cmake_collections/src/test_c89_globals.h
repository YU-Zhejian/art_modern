/* Does the compiler advertise C89 conformance?
   Do not test the value of __STDC__, because some compilers set it to 0
   while being otherwise adequately conformant. */
#if !defined __STDC__
#error "Compiler does not advertise C89 conformance"
#endif

#include <stdarg.h>
#include <stddef.h>
struct stat;
/* Most of the following tests are stolen from RCS 5.7 src/conf.sh.  */
struct buf {
    int x;
};
struct buf* (*rcsopen)(struct buf*, struct stat*, int);
static char* e(char** p, int i) { return p[i]; }
static char* f(char* (*g)(char**, int), char** p, ...)
{
    char* s;
    va_list v;
    va_start(v, p);
    s = g(p, va_arg(v, int));
    va_end(v);
    return s;
}

/* C89 style stringification. */
#define noexpand_stringify(a) #a
const char* stringified = noexpand_stringify(arbitrary + token = sequence);

/* C89 style token pasting.  Exercises some of the corner cases that
   e.g. old MSVC gets wrong, but not very hard. */
#define noexpand_concat(a, b) a##b
#define expand_concat(a, b) noexpand_concat(a, b)
/* extern int vA;*/
/* extern int vbee;*/
int vA;
int vbee;
#define aye A
#define bee B
int* pvA = &expand_concat(v, aye);
int* pvbee = &noexpand_concat(v, bee);

/* OSF 4.0 Compaq cc is some sort of almost-ANSI by default.  It has
   function prototypes and stuff, but not \xHH hex character constants.
   These do not provoke an error unfortunately, instead are silently treated
   as an "x".  The following induces an error, until -std is added to get
   proper ANSI mode.  Curiously \x00 != x always comes out true, for an
   array size at least.  It is necessary to write \x00 == 0 to get something
   that is true only with -std.  */
int osf4_cc_array['\x00' == 0 ? 1 : -1];

/* IBM C 6 for AIX is almost-ANSI by default, but it replaces macro parameters
   inside strings and character constants.  */
#define FOO(x) 'x'
int xlc6_cc_array[FOO(a) == 'x' ? 1 : -1];

int test(int i, double x);
struct s1 {
    int (*f)(int a);
};
struct s2 {
    int (*f)(double a);
};
int pairnames(int, char**, int* (*)(struct buf*, struct stat*, int), int, int);

/** Copied from GNU Autoconf a6a47a2f
 Modification: Return 0 and 1 reverted.
 */
int c89_const_main(void)
{
#ifndef __cplusplus
    /* Ultrix mips cc rejects this sort of thing.  */
    typedef int charset[2];
    const charset cs = { 0, 0 };
    /* SunOS 4.1.1 cc rejects this.  */
    char const* const* pcpcc;
    char** ppc;
    /* NEC SVR4.0.2 mips cc rejects this.  */
    struct point {
        int x, y;
    };
    static struct point const zero = { 0, 0 };
    /* IBM XL C 1.02.0.0 rejects this.
       It does not let you subtract one const X* pointer from another in
       an arm of an if-expression whose if-part is not a constant
       expression */
    const char* g = "string";
    pcpcc = &g + (g ? g - g : 0);
    /* HPUX 7.0 cc rejects these. */
    ++pcpcc;
    ppc = (char**)pcpcc;
    pcpcc = (char const* const*)ppc;
    { /* SCO 3.2v4 cc rejects this sort of thing.  */
        char tx;
        char* t = &tx;
        char const* s = 0 ? (char*)0 : (char const*)0;

        *t++ = 0;
        if (s)
            return 1;
    }
    { /* Someone thinks the Sun supposedly-ANSI compiler will reject this.  */
        int x[] = { 25, 17 };
        const int* foo = &x[0];
        ++foo;
    }
    { /* Sun SC1.0 ANSI compiler rejects this -- but not the above. */
        typedef const int* iptr;
        iptr p = 0;
        ++p;
    }
    { /* IBM XL C 1.02.0.0 rejects this sort of thing, saying
         "k.c", line 2.27: 1506-025 (S) Operand must be a modifiable lvalue. */
        struct s {
            int j;
            const int* ap[3];
        } bx;
        struct s* b = &bx;
        b->j = 5;
    }
    { /* ULTRIX-32 V3.1 (Rev 9) vcc rejects this */
        const int foo = 10;
        if (!foo)
            return 1;
    }
    return !(!cs[0] && !zero.x);
#endif
}

int c89main(int argc, char** argv)
{
    int ok = 0;
    ok |= (argc == 0 || f(e, argv, 0) != argv[0] || f(e, argv, 1) != argv[1]);
    ok |= c89_const_main();
    return ok;
}
