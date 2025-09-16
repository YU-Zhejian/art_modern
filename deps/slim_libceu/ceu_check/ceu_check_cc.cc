#include "ceu_check/ceu_check_cc.hh"
#include <sstream>

std::string ceu_check_get_compiler_info()
{
    std::ostringstream oss;
#if defined(__DATE__)
    std::string date_str = __DATE__;
#else
    std::string date_str = "unknown date";
#endif
#if defined(__TIME__)
    std::string time_str = __TIME__;
#else
    std::string time_str = "unknown time";
#endif
    oss << "Compiled at: " << date_str << ", " << time_str << std::endl;
    oss << "Compiler Identification:" << std::endl;
#if defined(CEU_COMPILER_IS_INTEL_CLANG)
    oss << "\t" << "Intel Clang compatible version number: " << __INTEL_CLANG_COMPILER / 10000 << '.'
        << __INTEL_CLANG_COMPILER % 10000 / 100 << '.' << __INTEL_CLANG_COMPILER % 100 << std::endl;
#endif
#if defined(CEU_COMPILER_IS_AOCC)
    oss << "\t" << "AMD Optimizing C++ Compiler (AOCC) compatible version number: "
#ifdef __aocc_major__
        << __aocc_major__ << '.' << __aocc_minor__ << '.' << __aocc_patchlevel__
#else
        << "unknown"
#endif
        << std::endl;
#endif
#if defined(CEU_COMPILER_IS_ICC)
    oss << "\t" << "Intel Compiler Classic (ICC) compatible version number: " << std::to_string(__ICC)

#ifdef __INTEL_COMPILER_UPDATE
        << "." << std::to_string(__INTEL_COMPILER_UPDATE)
#endif
        << std::endl;
#endif

#if defined(CEU_COMPILER_IS_MSVC)
#ifdef _MSC_BUILD
    std::string msc_build_ver = std::to_string(_MSC_BUILD);
#else
    std::string msc_build_ver = "unknown";
#endif

#ifdef _MSC_FULL_VER
    std::string msc_internal_ver = std::to_string(_MSC_FULL_VER % 100000);
#else
    std::string msc_internal_ver = "unknown";
#endif

#ifdef _MSC_VER
    std::string msc_major_ver = std::to_string(_MSC_VER / 100);
    std::string msc_minor_ver = std::to_string(_MSC_VER % 100);
#else
    std::string msc_major_ver = "unknown";
    std::string msc_minor_ver = "unknown";
#endif
    oss << "\t" << "MSVC compatible version number: " << msc_major_ver << '.' << msc_minor_ver << '.'
        << msc_internal_ver << '.' << msc_build_ver << ", with Visual Studio ver. " << CEU_VISUAL_STUDIO_VER << std::endl;
#endif

#if defined(CEU_COMPILER_IS_NVHPC)
    oss << "\t" << "NVidia High-Performance Compiler (NVHPC) compatible version number: " << __NVCOMPILER_MAJOR__ << '.'
        << __NVCOMPILER_MINOR__ << '.' << __NVCOMPILER_PATCHLEVEL__ << std::endl;
#endif
#if defined(CEU_COMPILER_IS_PGIC)
    oss << "\t" << "PGI Compiler (PGIC) compatible version number: " << __PGIC__ << '.' << __PGIC_MINOR__ << '.'
        << __PGIC_PATCHLEVEL__ << std::endl;
#endif

#if defined(CEU_COMPILER_IS_TCC)
    oss << "\t" << "Tiny C Compiler (TCC) compatible version number: ";
    int tcc_major = __TINYC__ / 10000;
    int tcc_minor = (__TINYC__ - tcc_major * 10000) / 100;
    int tcc_patchlevel = __TINYC__ % 100;
    oss << tcc_major << '.' << tcc_minor << '.' << tcc_patchlevel << std::endl;
#endif
#if defined(CEU_COMPILER_IS_EDG)
    int edg_major = __EDG_VERSION__ / 100;
    int edg_minor = (__EDG_VERSION__ - edg_major * 100);
    oss << "\t" << "EDG compatible version number: " << edg_major << '.' << edg_minor << std::endl;
#endif

#if defined(CEU_COMPILER_IS_BORLAND)
    int borland_major = __BORLANDC__ / 256;
    int borland_revision = (__BORLANDC__ - 256 * borland_major) / 16 + __BORLANDC__ % 16;
    oss << "\t" << "Borland compatible version number: " << borland_major << '.' << borland_revision << ", with "
        << CEU_CPPB_VERSION << std::endl;
#endif
#if defined(CEU_COMPILER_IS_CLANG)
    oss << "\t" << "Clang compatible version number: ";
#ifdef __clang_major__
    oss << __clang_major__ << '.' << __clang_minor__ << '.' << __clang_patchlevel__;
#else
    oss << "UNKNOWN";
#endif
#ifdef __clang_version__
    oss << " (" << __clang_version__ << ")";
#endif
    oss << std::endl;
#endif

#if defined(CEU_COMPILER_IS_GCC)
    oss << "\t" << "GCC compatible version number: ";
#ifdef __GNUC_PATCHLEVEL__
    oss << __GNUC__ << '.' << __GNUC_MINOR__ << '.' << __GNUC_PATCHLEVEL__;
#else
    oss << __GNUC__ << '.' << __GNUC_MINOR__;
#endif
    oss << std::endl;
#endif
#if defined(__VERSION__)
    oss << "\t" << "__VERSION__ string: " << __VERSION__ << std::endl;
#endif
    return oss.str();
}
