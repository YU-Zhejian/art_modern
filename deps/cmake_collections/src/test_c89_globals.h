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
static char* e(char** p, int i)
{
    return p[i];
}
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
// extern int vA;
// extern int vbee;
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
int pairnames(int, char**, int* (*)(struct buf*, struct stat*, int),
    int, int);

int c89main(int argc, char** argv)
{
    int ok;
    ok |= (argc == 0 || f(e, argv, 0) != argv[0] || f(e, argv, 1) != argv[1]);
    return ok;
}
