#include "ceu_check/ceu_check_cxx_stdlib.hh"

#include "ceu_check/ceu_check_cxx_stdlib_macro.hh"
#include "ceu_check/ceu_constants.h"

#include <ostream>
#include <sstream>
#include <string>

std::string ceu_interpret_cxx_stdlib()
{
    std::ostringstream oss;
    oss << "C++ Standard Library: " << CEU_CXX_STDLIB_NAME << std::endl;
#ifdef CEU_CXX_STDLIB_IS_LIBCXX
    /** Determine whether the libcxx is before 16.0.0 */
    if (_LIBCPP_VERSION < 160000) {
        // The libcxx version is XX.Y.ZZ
        oss << "\tVersion: " << _LIBCPP_VERSION / 1000 << "." << (_LIBCPP_VERSION / 100) % 10 << "."
            << _LIBCPP_VERSION % 100 << " (pre-16 spec)" << std::endl;
    } else {
        // The libcxx version is XX.YY.ZZ
        oss << "\tVersion: " << _LIBCPP_VERSION / 10000 << "." << (_LIBCPP_VERSION / 100) % 100 << "."
            << _LIBCPP_VERSION % 100 << " (post-16 spec)" << std::endl;
    }
#elif defined(CEU_CXX_STDLIB_IS_LIBSTDCXX)
    oss << "\tVersion: " << CEU_CXX_STDLIB_GLIBCXX_VERSION << std::endl;
    oss << "\t_GLIBCXX_RELEASE: " <<
#ifdef _GLIBCXX_RELEASE
        _GLIBCXX_RELEASE
#else
        CEU_UNDEFINED
#endif
        << std::endl;
    oss << "\t_GLIBCXX_VERSION/_GLIBCPP_VERSION: " <<
#if defined(_GLIBCPP_VERSION) || defined(_GLIBCXX_VERSION)
#ifdef _GLIBCXX_VERSION
        _GLIBCXX_VERSION
#else
        _GLIBCPP_VERSION
#endif
#else
        CEU_UNDEFINED
#endif
        << std::endl;
#elif defined(CEU_CXX_STDLIB_IS_MSVC)
    oss << "\tVersion: " << _CPPLIB_VER / 100 << "." << _CPPLIB_VER % 100 << std::endl;
    oss << "\t\t_MSVC_STL_VERSION=" <<
#ifdef _MSVC_STL_VERSION
        _MSVC_STL_VERSION
#else
        CEU_UNDEFINED
#endif
        << std::endl;
    oss << "\t\t_MSVC_STL_UPDATE=" <<
#ifdef _MSVC_STL_UPDATE
        _MSVC_STL_UPDATE
#else
        CEU_UNDEFINED
#endif
        << std::endl;
#endif
    return oss.str();
}
