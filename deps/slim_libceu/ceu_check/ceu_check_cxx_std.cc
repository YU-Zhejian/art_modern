#include "ceu_check/ceu_check_cxx_std.hh"

#include "ceu_check/ceu_check_c_cxx_std_macro.h"

#include <sstream>
#include <string>

std::string ceu_interpret_cxx_std_version()
{
    std::string cxxstd_macro;
    std::ostringstream oss;
    oss << "Compile-time C++ std.: ver. " << CEU_CXX_STD;
#ifdef CEU_CXX_STD_VERSION_MACRO
    cxxstd_macro = std::to_string(CEU_CXX_STD_VERSION_MACRO);
#else
    cxxstd_macro = "_MSVC_LANG and __cplusplus undefined";
#endif
    oss << " (" << cxxstd_macro << ")" << std::endl;
    return oss.str();
}
