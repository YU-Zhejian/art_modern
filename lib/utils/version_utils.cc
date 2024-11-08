#include "version_utils.hh"
#include "art_modern_config.h"
#include "art_modern_constants.hh"

#include <iostream>

// Boost
#include <boost/algorithm/string/trim.hpp>
#include <boost/version.hpp>

// HTSLib
#include <htslib/hfile.h>
#include <htslib/hts.h>

// GSL
#ifdef USE_GSL_RANDOM
#include <gsl/gsl_version.h>
#endif

// MKL
#ifdef USE_ONEMKL_RANDOM
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

namespace labw::art_modern {
void print_htslib_version()
{
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
}

void print_boost_version()
{
    int patch_level = BOOST_VERSION % 100;
    int minor = BOOST_VERSION / 100 % 1000;
    int major = BOOST_VERSION / 100000;
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

void print_onemkl_version()
{
#ifdef USE_ONEMKL_RANDOM
    std::cout << "MKL Version: " << __INTEL_MKL__ << "." << __INTEL_MKL_MINOR__ << "." << __INTEL_MKL_UPDATE__ << "."
              << __INTEL_MKL_PATCH__ << " (" << INTEL_MKL_VERSION << ")" << std::endl;
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

void print_version()
{
    std::cout << "ART: " << ART_VERSION << ", ART_MODERN: " << ART_MODERN_VERSION << std::endl;
    std::cout << "ART_MODERN_LINK_LIBS: " << ART_MODERN_LINK_LIBS << std::endl;
    print_htslib_version();
    print_boost_version();
    print_gsl_version();
    print_onemkl_version();
    print_mpi_version();
    print_protobuf_version();
}
} // art_modern
// labw