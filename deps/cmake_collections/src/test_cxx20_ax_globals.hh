#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 202002L && !defined _MSC_VER

#error "This is not a C++20 compiler"

#else

#include <version>

namespace cxx20 {

// As C++20 supports feature test macros in the standard, there is no
// immediate need to actually test for feature availability on the
// Autoconf side.

} // namespace cxx20

#endif // __cplusplus < 202002L && !defined _MSC_VER
