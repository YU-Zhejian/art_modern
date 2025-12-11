#pragma once
#include <cstdlib> // NOLINT
#if defined(_LIBCPP_VERSION)
// For LLVM libc++
#define CEU_CXX_STDLIB_IS_LIBCXX
#define CEU_CXX_STDLIB_NAME "LLVM libc++"
#elif defined(__GLIBCPP__) || defined(__GLIBCXX__)
#ifdef __GLIBCXX__
#define CEU_CXX_STDLIB_GLIBCXX_VERSION __GLIBCXX__
#else
#define CEU_CXX_STDLIB_GLIBCXX_VERSION __GLIBCPP__
#endif
#define CEU_CXX_STDLIB_IS_LIBSTDCXX
#define CEU_CXX_STDLIB_NAME "GNU libstdc++"
#elif defined(_CPPLIB_VER)
#define CEU_CXX_STDLIB_IS_MSVC
#define CEU_CXX_STDLIB_NAME "MSVC STL (by P. J. Plauger/Dinkumware, Ltd.)"
#else
#define CEU_CXX_STDLIB_IS_UNKNOWN
#define CEU_CXX_STDLIB_NAME "Unknown C++ Standard Library"
#endif
