#include "libceu_stddef.h" /* NOLINT */

#include "ceu_check/ceu_check_cxx_compiler.hh"

#include "ceu_check/ceu_check_c_cxx_compiler_macro.h"
#include "ceu_utils/ceu_constants.h" /* NOLINT: for CEU_UNDEFINED */

#include "ceu_utils/ceu_compiler_names.h"

#include <ostream>
#include <sstream>
#include <string>

std::string ceu_check_get_cxx_compiler_info()
{
    std::ostringstream oss;
#if defined(CEU_REPRODUCIBLE_BUILDS) && (CEU_REPRODUCIBLE_BUILDS == 1)
    oss << "C++ Compiled at: N/A due to reproducible build" << std::endl;
#else
    const std::string date_str =
#if defined(__DATE__)
        __DATE__;
#else
        CEU_UNDEFINED;
#endif
    const std::string time_str =
#if defined(__TIME__)
        __TIME__;
#else
        CEU_UNDEFINED;
#endif
    oss << "C++ Compiled at: " << date_str << ", " << time_str << std::endl;
#endif
    oss << "C++ Compiler Identification: " << CEU_COMPILER_NAME << std::endl;

#if defined(CEU_COMPILER_IS_TI)
    long ti_version =
#ifdef __TI_COMPILER_VERSION__
        __TI_COMPILER_VERSION__
#elif defined(__TI_COMPILER_VERSION)
        __TI_COMPILER_VERSION
#else
        0L
#endif
        ;
    oss << "\t" << CEU_COMPILER_NAME_TI << " compatible version number: " << ti_version / 1000000L << '.'
        << (ti_version % 1000000L) / 1000 << '.' << ti_version % 1000 << std::endl;
#endif

#if defined(CEU_COMPILER_IS_INTEL_CLANG)
    oss << "\t" << CEU_COMPILER_NAME_INTEL_ONEAPI_DPCPP
        << " compatible version number: " << __INTEL_CLANG_COMPILER / 10000 << '.'
        << __INTEL_CLANG_COMPILER % 10000 / 100 << '.' << __INTEL_CLANG_COMPILER % 100 << std::endl;
#endif
#if defined(CEU_COMPILER_IS_ARM_COMPILER_LINUX)
    oss << "\t" << CEU_COMPILER_NAME_ARM_COMPILER_LINUX << " compatible version number: " << __armclang_major__ << '.'
        << __armclang_minor__ << " (" << __armclang_version__ << ")" << " build " << __ARM_LINUX_COMPILER_BUILD__
        << std::endl;
#endif
#if defined(CEU_COMPILER_IS_ARM_COMPILER_EMBEDDED)
    oss << "\t" << CEU_COMPILER_NAME_ARM_COMPILER_EMBEDDED << " compatible version number: " <<
#ifdef __ARMCC_VERSION
        // P.VV.BBBB
        __ARMCC_VERSION / 1000000 << "." << (__ARMCC_VERSION / 10000) % 100 << "." << __ARMCC_VERSION % 10000
#else
        "UNDEFINED"
#endif
        << std::endl;
#endif

#if defined(CEU_COMPILER_IS_AOCC)
    oss << "\t" << CEU_COMPILER_NAME_AOCC << " compatible version number: "
#ifdef __aocc_major__
        << __aocc_major__ << '.' << __aocc_minor__ << '.' << __aocc_patchlevel__
#else
        << CEU_UNDEFINED
#endif
        << std::endl;
#endif
#if defined(CEU_COMPILER_IS_ICC)
    oss << "\t" << CEU_COMPILER_NAME_ICC << " compatible version number: " << std::to_string(__ICC)

#ifdef __INTEL_COMPILER_UPDATE
        << "." << std::to_string(__INTEL_COMPILER_UPDATE)
#endif
        << std::endl;
#endif

#if defined(CEU_COMPILER_IS_MSVC)
    const std::string msc_build_ver =
#ifdef _MSC_BUILD
        std::to_string(_MSC_BUILD);
#else
        CEU_UNDEFINED;
#endif

    const std::string msc_internal_ver =
#ifdef _MSC_FULL_VER
        std::to_string(_MSC_FULL_VER % 100000);
#else
        CEU_UNDEFINED;
#endif

#ifdef _MSC_VER
    const std::string msc_major_ver = std::to_string(_MSC_VER / 100);
    const std::string msc_minor_ver = std::to_string(_MSC_VER % 100);
#else
    const std::string msc_major_ver = CEU_UNDEFINED;
    const std::string msc_minor_ver = CEU_UNDEFINED;
#endif
    oss << "\t" << CEU_COMPILER_NAME_MSVC << " compatible version number: " << msc_major_ver << '.' << msc_minor_ver
        << '.' << msc_internal_ver << '.' << msc_build_ver << ", with Visual Studio ver. " << CEU_VISUAL_STUDIO_VER
        << std::endl;
#endif

#if defined(CEU_COMPILER_IS_NVHPC)
    oss << "\t" << CEU_COMPILER_NAME_NVHPC << " compatible version number: " << __NVCOMPILER_MAJOR__ << '.'
        << __NVCOMPILER_MINOR__ << '.' << __NVCOMPILER_PATCHLEVEL__ << std::endl;
#endif
#if defined(CEU_COMPILER_IS_PGIC)
    oss << "\t" << CEU_COMPILER_NAME_PGIC << " compatible version number: " << __PGIC__ << '.' << __PGIC_MINOR__ << '.'
        << __PGIC_PATCHLEVEL__ << std::endl;
#endif

#if defined(CEU_COMPILER_IS_TCC)
    oss << "\t" << CEU_COMPILER_NAME_TCC << " compatible version number: ";
    int tcc_major = __TINYC__ / 10000;
    int tcc_minor = (__TINYC__ - tcc_major * 10000) / 100;
    int tcc_patchlevel = __TINYC__ % 100;
    oss << tcc_major << '.' << tcc_minor << '.' << tcc_patchlevel << std::endl;
#endif
#if defined(CEU_COMPILER_IS_EDG)
    int edg_major = __EDG_VERSION__ / 100;
    int edg_minor = (__EDG_VERSION__ - edg_major * 100);
    oss << "\t" << CEU_COMPILER_NAME_EDG << " compatible version number: " << edg_major << '.' << edg_minor
        << std::endl;
#endif

#if defined(CEU_COMPILER_IS_BORLAND)
    int borland_major = __BORLANDC__ / 256;
    int borland_revision = (__BORLANDC__ - 256 * borland_major) / 16 + __BORLANDC__ % 16;
    oss << "\t" << CEU_COMPILER_NAME_BORLAND << "compatible version number: " << borland_major << '.'
        << borland_revision << ", with " << CEU_CPPB_VERSION << std::endl;
    oss << "\t\t__CODEGEARC_VERSION__="
#ifdef __CODEGEARC_VERSION__
        << std::hex << ((__CODEGEARC_VERSION__ & 0xFF000000) >> 24) << '.'
        << ((__CODEGEARC_VERSION__ & 0x00FF0000) >> 16) << '.' << std::dec << ((__CODEGEARC_VERSION__ & 0x0000FFFF))
        << std::endl;
#else
        << "=UNDEFINED" << std::endl;
#endif
#endif
#if defined(CEU_COMPILER_IS_CLANG)
    oss << "\t" << CEU_COMPILER_NAME_CLANG << " compatible version number: ";
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
    oss << "\t" << CEU_COMPILER_NAME_GCC << " compatible version number: ";
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
