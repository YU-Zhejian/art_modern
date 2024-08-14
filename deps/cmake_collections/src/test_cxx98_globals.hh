// Does the compiler advertise C++98 conformance?
#if !defined __cplusplus || __cplusplus < 199711L
#error "Compiler does not advertise C++98 conformance"
#endif

// These inclusions are to reject old compilers that
// lack the unsuffixed header files.
#include <cstdlib>
#include <exception>

// FIXME: Following are the original GNU Autoconf files
// <cassert> and <cstring> are *not* freestanding headers in C++98.
// extern void assert (int);
// namespace std {
//     extern int strcmp (const char *, const char *);
// }
#include <cassert>
#include <cstring>

// Namespaces, exceptions, and templates were all added after "C++ 2.0".
using std::exception;
using std::strcmp;

namespace {

void test_exception_syntax()
{
    try {
        throw "test";
    } catch (const char* s) {
        // Extra parentheses suppress a warning when building autoconf itself,
        // due to lint rules shared with more typical C programs.
        assert(!(strcmp)(s, "test"));
    }
}

template <typename T>
struct test_template {
    T const val;
    explicit test_template(T t)
        : val(t)
    {
    }
    template <typename U>
    T add(U u) { return static_cast<T>(u) + val; }
};

} // anonymous namespace

int cxx98main(int argc, char** argv)
{
    int ok = 0;
    assert(argc);
    // assert (! argv[0]);
    {
        test_exception_syntax();
        test_template<double> tt(2.0);
        assert(tt.add(4) == 6.0);
        assert(true && !false);
    }
    return ok;
}
