#include <iostream>

class Myclass {
public:
    Myclass(Myclass&&) noexcept { std::cout << "Move constructor called" << std::endl; };
    Myclass& operator=(Myclass&&) = delete;
    Myclass(const Myclass&) noexcept { std::cout << "Copy constructor called" << std::endl; };
    Myclass& operator=(const Myclass&) = delete;

    Myclass() { std::cout << "Myclass()" << std::endl; }
    ~Myclass() { std::cout << "~Myclass()" << std::endl; }
    void foo() const { std::cout << "foo()" << std::endl; }
};

Myclass f() { return Myclass(); }

int main(int, char**)
{
    Myclass a = f();
    a.foo();
    return 0;
}
