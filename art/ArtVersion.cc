#include <art_modern_config.h>
#include <boost/version.hpp>
#include <htslib/hts.h>
#include <iostream>

#include "ArtVersion.hh"
#include "art_modern_constants.hh"
#include "htslib/hfile.h"

#ifdef USE_GSL_RANDOM
#include <gsl/gsl_version.h>
#endif

#ifdef USE_ONEMKL_RANDOM
#include <mkl_version.h>
#endif

#include <mpi.h>

namespace labw {
namespace art_modern {
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
        std::cout << "MKL Version: " << __INTEL_MKL__ << "." << __INTEL_MKL_MINOR__ << "." << __INTEL_MKL_UPDATE__
                  << "." << __INTEL_MKL_PATCH__ << " (" << INTEL_MKL_VERSION << ")" << std::endl;
#else
        std::cout << "MKL: not used" << std::endl;
#endif
    }

    void print_mpi_version()
    {
#ifdef WITH_MPI
        int major;
        int minor;
        MPI_Get_version(&major, &minor);

        char libversion[MPI_MAX_LIBRARY_VERSION_STRING];
        int resultlen;
        MPI_Get_library_version(libversion, &resultlen);
        std::cout << "MPI:" << std::endl;
        std::cout << "\tStandard Version: " << major << "." << minor << std::endl;
        std::cout << "\tLibrary Version: " << std::string(libversion, resultlen) << std::endl;

        char vendor_string[MPI_MAX_PROCESSOR_NAME];
        int vendor_len;
        MPI_Get_processor_name(vendor_string, &vendor_len);

        std::cout << "\tVendor: " << std::string(vendor_string, vendor_len) << std::endl;
#else
        std::cout << "MPI: not used" << std::endl;
#endif
    }

    void print_version()
    {
        std::cout << "ART: " << ART_VERSION << ", ART_MODERN: " << ART_MODERN_VERSION << std::endl;
        print_htslib_version();
        print_boost_version();
        print_gsl_version();
        print_onemkl_version();
        print_mpi_version();
    }
} // art_modern
} // labw