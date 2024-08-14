#include "test_cxx98_globals.hh"
// Does the compiler advertise C++ 2011 conformance?
#if !defined __cplusplus || __cplusplus < 201103L
#error "Compiler does not advertise C++11 conformance"
#endif

namespace cxx11test {
constexpr int get_val() { return 20; }

struct testinit {
    int i;
    double d;
};

class delegate {
public:
    delegate(int n)
        : n(n)
    {
    }
    delegate()
        : delegate(2354)
    {
    }

    virtual int getval() { return this->n; };

protected:
    int n;
};

class overridden : public delegate {
public:
    overridden(int n)
        : delegate(n)
    {
    }
    virtual int getval() override final { return this->n * 2; }
};

class nocopy {
public:
    nocopy(int i)
        : i(i)
    {
    }
    nocopy() = default;
    nocopy(const nocopy&) = delete;
    nocopy& operator=(const nocopy&) = delete;

private:
    int i;
};

// for testing lambda expressions
template <typename Ret, typename Fn>
Ret eval(Fn f, Ret v)
{
    return f(v);
}

// for testing variadic templates and trailing return types
template <typename V>
auto sum(V first) -> V
{
    return first;
}
template <typename V, typename... Args>
auto sum(V first, Args... rest) -> V
{
    return first + sum(rest...);
}
}
int cxx11main(int argc, char** argv)
{
    int ok = 0;
    {
        // Test auto and decltype
        auto a1 = 6538;
        auto a2 = 48573953.4;
        auto a3 = "String literal";

        int total = 0;
        for (auto i = a3; *i; ++i) {
            total += *i;
        }

        decltype(a2) a4 = 34895.034;
    }
    {
        // Test constexpr
        short sa[cxx11test::get_val()] = { 0 };
    }
    {
        // Test initializer lists
        cxx11test::testinit il = { 4323, 435234.23544 };
    }
    {
        // Test range-based for
        int array[] = { 9, 7, 13, 15, 4, 18, 12, 10, 5, 3,
            14, 19, 17, 8, 6, 20, 16, 2, 11, 1 };
        for (auto& x : array) {
            x += 23;
        }
    }
    {
        // Test lambda expressions
        using cxx11test::eval;
        assert(eval([](int x) { return x * 2; }, 21) == 42);
        double d = 2.0;
        assert(eval([&](double x) { return d += x; }, 3.0) == 5.0);
        assert(d == 5.0);
        assert(eval([=](double x) mutable { return d += x; }, 4.0) == 9.0);
        assert(d == 5.0);
    }
    {
        // Test use of variadic templates
        using cxx11test::sum;
        auto a = sum(1);
        auto b = sum(1, 2);
        auto c = sum(1.0, 2.0, 3.0);
    }
    {
        // Test constructor delegation
        cxx11test::delegate d1;
        cxx11test::delegate d2();
        cxx11test::delegate d3(45);
    }
    {
        // Test override and final
        cxx11test::overridden o1(55464);
    }
    {
        // Test nullptr
        char* c = nullptr;
    }
    {
        // Test template brackets
        test_template<::test_template<int>> v(test_template<int>(12));
    }
    {
        // Unicode literals
        char const* utf8 = u8"UTF-8 string \u2500";
        char16_t const* utf16 = u"UTF-8 string \u2500";
        char32_t const* utf32 = U"UTF-32 string \u2500";
    }
    return ok;
}