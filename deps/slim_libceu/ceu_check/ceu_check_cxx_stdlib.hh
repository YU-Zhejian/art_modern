#pragma once
#include <cstdlib> // NOLINT: for _LIBCPP_VERSION/__GLIBCXX__/_MSC_VER
// See also: https://en.cppreference.com/w/cpp/header/ciso646.html
// DO NOT USE <ciso646> or <iso646.h> to detect standard library!
// See also: <https://github.com/microsoft/STL/blob/main/stl/inc/yvals_core.h>
// L936-938
// #define _CPPLIB_VER       650
// #define _MSVC_STL_VERSION 145
// #define _MSVC_STL_UPDATE  202511L
// See also: https://github.com/cpredef/predef/blob/master/Libraries.md

#if defined(_LIBCPP_VERSION)
// For LLVM libc++
#define CEU_CXX_STDLIB_IS_LIBCXX
#elif defined(__GLIBCPP__) || defined(__GLIBCXX__)
#define CEU_CXX_STDLIB_IS_LIBSTDCXX
#elif defined(_CPPLIB_VER)
#define CEU_CXX_STDLIB_IS_MSVC
#else
#define CEU_CXX_STDLIB_IS_UNKNOWN
#endif
