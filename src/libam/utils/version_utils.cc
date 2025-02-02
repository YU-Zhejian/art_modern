#include "art_modern_config.h"

#include "libam/utils/version_utils.hh"

#include "libam/Constants.hh"

#include "ceu_check/ceu_check_c_cxx_std.hh"
#include "ceu_check/ceu_check_cc.hh"
#include "ceu_check/ceu_check_ctypes_limit.hh"
#include "ceu_check/ceu_check_os.hh"

// Boost
#include <boost/algorithm/string/join.hpp>
#include <boost/version.hpp>

// HTSLib
#include <htslib/hfile.h>
#include <htslib/hts.h>

#ifdef USE_BS_PARALLEL
#include <BS_thread_pool.hpp>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

// GSL
#ifdef WITH_GSL
#include <gsl/gsl_version.h>
#endif

// MKL
#ifdef WITH_ONEMKL
#include <mkl_version.h>
#endif

// MPI
#ifdef WITH_MPI
#include <mpi.h>
#endif

// Protobuf
#ifdef WITH_PROTOBUF
#include <google/protobuf/stubs/common.h>
#endif

// malloc
#if defined(WITH_MIMALLOC)
#include <mimalloc.h>
#elif defined(WITH_JEMALLOC)
#include <jemalloc/jemalloc.h>
#endif

#include <fmt/core.h>

// CPPSTDLIB
#include <iostream>
#include <string>
#include <vector>

namespace labw::art_modern {
namespace {
    void print_htslib_version()
    {
        // NOLINTBEGIN
        std::cout << "USING HTSLib: " << USED_HTSLIB_NAME << ", ver. " << hts_version() << std::endl;
        std::cout << "\tFeatures: " << hts_feature_string() << std::endl;
        std::cout << "\tCC: " << hts_test_feature(HTS_FEATURE_CC) << std::endl;
        std::cout << "\tCPPFLAGS: " << hts_test_feature(HTS_FEATURE_CPPFLAGS) << std::endl;
        std::cout << "\tCFLAGS: " << hts_test_feature(HTS_FEATURE_CFLAGS) << std::endl;
        std::cout << "\tLDFLAGS: " << hts_test_feature(HTS_FEATURE_LDFLAGS) << std::endl;
        std::cout << "\tHTSlib URL scheme handlers present:" << std::endl;

        // Following code were modified from SAMTtools ver. 1.13.
        const char* plugins[100];
        int np = 100;

        if (hfile_list_plugins(plugins, &np) >= 0) {
            for (int i = 0; i < np; i++) {
                const char* sc_list[100];
                int nschemes = 100;
                if (hfile_list_schemes(plugins[i], sc_list, &nschemes) < 0) {
                    continue;
                }
                std::cout << "\t\t" << plugins[i] << ": ";
                if (nschemes == 0) {
                    std::cout << std::endl;
                } else {
                    for (int j = 0; j < nschemes - 1; j++) {
                        std::cout << sc_list[j] << ", ";
                    }
                }
                std::cout << sc_list[nschemes - 1];
                std::cout << std::endl;
            }
        }
        // NOLINTEND
    }

    void print_boost_version()
    {
        constexpr int patch_level = BOOST_VERSION % 100;
        constexpr int minor = BOOST_VERSION / 100 % 1000;
        constexpr int major = BOOST_VERSION / 100000;
        std::cout << "BOOST: " << major << "." << minor << "." << patch_level << std::endl;
    }

    void print_gsl_version()
    {
#ifdef USE_GSL_RANDOM
        std::cout << "GSL: " << gsl_version << std::endl;
#else
        std::cout << "GSL: not used" << std::endl;
#endif
    }

    void print_bs_version()
    {
#ifdef USE_BS_PARALLEL
        std::cout << "BS::thread_pool: " << BS_THREAD_POOL_VERSION_MAJOR << "." << BS_THREAD_POOL_VERSION_MINOR << "."
                  << BS_THREAD_POOL_VERSION_PATCH << std::endl;
#else
        std::cout << "BS::thread_pool: not used" << std::endl;
#endif
    }

    void print_onemkl_version()
    {
#ifdef USE_ONEMKL_RANDOM
        std::cout << "MKL Version: " << __INTEL_MKL__ << "." << __INTEL_MKL_MINOR__ << "." << __INTEL_MKL_UPDATE__
                  << "." << __INTEL_MKL_PATCH__ << " (" << INTEL_MKL_VERSION << ")" << std::endl;
#else
        std::cout << "MKL: not used" << std::endl;
#endif
    }

    void print_mpi_version()
    {
#ifdef WITH_MPI
        std::cout << "MPI:" << std::endl;

        int major;
        int minor;
        MPI_Get_version(&major, &minor);
        std::cout << "\tStandard Version: " << major << "." << minor << std::endl;

        char libversion_string[MPI_MAX_LIBRARY_VERSION_STRING];
        int libversion_len;
        MPI_Get_library_version(libversion_string, &libversion_len);
        auto lib_ver_str = std::string(libversion_string, libversion_len);
        boost::algorithm::trim(lib_ver_str);
        std::cout << "\tLibrary Version: " << lib_ver_str << std::endl;
#else
        std::cout << "MPI: not used" << std::endl;
#endif
    }

    void print_protobuf_version()
    {
#ifdef WITH_PROTOBUF
        const int version = GOOGLE_PROTOBUF_VERSION;

        // Extract major, minor, and patch numbers
        const int major = version / 1000000;
        const int minor = (version % 1000000) / 1000;
        const int patch = version % 1000;

        // Print the version
        std::cout << "Protobuf: " << major << "." << minor << "." << patch << std::endl;
#else
        std::cout << "Protobuf: not used" << std::endl;
#endif
    }

    void print_openmp_version()
    {
#ifdef WITH_OPENMP
        std::cout << "OpenMP (C): " << OpenMP_C_SPEC_DATE <<
#ifdef OpenMP_C_VERSION
            " v" << OpenMP_C_VERSION <<
#endif
            std::endl;
        std::cout << "OpenMP (CXX): " << OpenMP_CXX_SPEC_DATE <<
#ifdef OpenMP_CXX_VERSION
            " v" << OpenMP_CXX_VERSION <<
#endif
            std::endl;
        std::cout << "OpenMP Macros: _OPENMP=";
#ifdef _OPENMP
        try {
            std::cout << _OPENMP;
        } catch (std::exception& e) {
            std::cout << "YES";
        }
#else
        std::cout << "NO";
#endif
        std::cout << " _OPENMP_SIMD=";
#ifdef _OPENMP_SIMD
        try {
            std::cout << _OPENMP_SIMD;
        } catch (std::exception& e) {
            std::cout << "YES";
        }
        std::cout << "YES";
#else
        std::cout << "NO";
#endif
        std::cout << std::endl;

#else
        std::cout << "OpenMP: not used" << std::endl;
#endif
    }

    void print_simde_version()
    {
        std::vector<std::string> simd_info;
#if (0)
        std::cout << "SIMDE: " << SIMDE_VERSION_MAJOR << "." << SIMDE_VERSION_MINOR << "." << SIMDE_VERSION_MICRO
                  << std::endl;
#ifdef SIMDE_X86_SSE_NATIVE
        simd_info.emplace_back("SSE");
#endif
#ifdef SIMDE_X86_SSE2_NATIVE
        simd_info.emplace_back("SSE2");
#endif
        std::cout << "\twith ISE: " << boost::algorithm::join(simd_info, " ") << std::endl;
#else
        std::cout << "SIMDE: N/A" << std::endl;
#ifdef __MMX__
        simd_info.emplace_back("MMX");
#endif
#ifdef __SSE2__
        simd_info.emplace_back("SSE2");
#endif
#ifdef __AVX2__
        simd_info.emplace_back("AVX2");
#endif
        std::cout << "\twith ISE: " << boost::algorithm::join(simd_info, " ") << std::endl;
#endif
    }

    void print_malloc_version()
    {
#if defined(WITH_MIMALLOC)
        std::cout << "mimalloc: " << MI_MALLOC_VERSION / 100 << "." << MI_MALLOC_VERSION % 100 << std::endl;
#elif defined(WITH_JEMALLOC)
#ifdef JEMALLOC_VERSION
        std::cout << "jemalloc: " << JEMALLOC_VERSION << std::endl;
#else
        std::cout << "jemalloc: " << JEMALLOC_VERSION_MAJOR << "." << JEMALLOC_VERSION_MINOR << "."
                  << JEMALLOC_VERSION_BUGFIX << std::endl;
#endif
#else
        std::cout << "*malloc: not used" << std::endl;
#endif
    }

    void print_fmt_version()
    {
        std::cout << "{fmt}: " << FMT_VERSION / 10000 << "." << FMT_VERSION / 100 % 100 << "." << FMT_VERSION % 100
                  << std::endl;
    }

} // namespace

void print_version()
{
    std::cout << "ART: " << ART_VERSION << ", ART_MODERN: " << ART_MODERN_VERSION << std::endl;
#ifdef WITH_GIT
    std::cout << "On git commit: (" << CEU_CM_GIT_COMMIT_HASH << ") "
              << CEU_CM_GIT_COMMIT_DATE << std::endl;
#else
    std::cout << "Git: N/A" << std::endl;
#endif
    std::cout << "ART_MODERN_LINK_LIBS: " << ART_MODERN_LINK_LIBS << std::endl;
    print_htslib_version();
    print_fmt_version();
    print_boost_version();
    print_gsl_version();
    print_onemkl_version();
    print_mpi_version();
    print_protobuf_version();
    print_openmp_version();
    print_simde_version();
    print_bs_version();
    print_malloc_version();
    std::cout << ceu_interpret_c_std_version();
    std::cout << ceu_interpret_cxx_std_version();
    std::cout << ceu_check_get_compiler_info();
    std::cout << ceu_check_get_ctypes_limit_info();
    std::cout << ceu_check_get_compile_time_os_info();
    std::cout << ceu_check_get_run_time_os_info();
}
} // namespace labw::art_modern
// labw