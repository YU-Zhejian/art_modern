#include "ceu_check/ceu_check_c_cxx_std.hh"
#include <sstream>

std::string ceu_interpret_c_std_version(void)
{
    std::ostringstream oss;
    std::string cstd_macro;
    oss << "Compile-time C std.: ver. " << CEU_C_STD;
#ifdef CEU_C_STD_VERSION_MACRO
    cstd_macro = std::to_string(CEU_C_STD_VERSION_MACRO);
#else
    cstd_macro = "__STDC_VERSION__ undefined";
#endif
    oss << " (" << cstd_macro << ")" << std::endl;
    return oss.str();
}

std::string ceu_interpret_cxx_std_version(void)
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
