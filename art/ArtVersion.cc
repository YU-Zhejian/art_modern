#include <art_modern_config.h>
#include <boost/version.hpp>
#include <htslib/hts.h>
#include <iostream>

#include "ArtVersion.hh"
#include "art_modern_constants.hh"
#include "htslib/hfile.h"

using namespace std;

namespace labw {
namespace art_modern {
    void print_htslib_version()
    {
        std::cout << "USING HTSLib: " << USED_HTSLIB_NAME << ", ver. " << hts_version() << endl;
        std::cout << "\tFeatures: " << hts_feature_string() << endl;
        std::cout << "\tCC: " << hts_test_feature(HTS_FEATURE_CC) << endl;
        std::cout << "\tCPPFLAGS: " << hts_test_feature(HTS_FEATURE_CPPFLAGS) << endl;
        std::cout << "\tCFLAGS: " << hts_test_feature(HTS_FEATURE_CFLAGS) << endl;
        std::cout << "\tLDFLAGS: " << hts_test_feature(HTS_FEATURE_LDFLAGS) << endl;
        std::cout << "\tHTSlib URL scheme handlers present:" << endl;

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
                    std::cout << endl;
                } else {
                    for (int j = 0; j < nschemes - 1; j++) {
                        std::cout << sc_list[j] << ", ";
                    }
                }
                std::cout << sc_list[nschemes - 1];
                std::cout << endl;
            }
        }
    }

    void print_boost_version()
    {
        int patch_level = BOOST_VERSION % 100;
        int minor = BOOST_VERSION / 100 % 1000;
        int major = BOOST_VERSION / 100000;
        std::cout << "BOOST: " << major << "." << minor << "." << patch_level << endl;
    }

    void print_version()
    {
        std::cout << "ART: " << ART_VERSION << ", ART_MODERN: " << ART_MODERN_VERSION << endl;
        print_htslib_version();
        print_boost_version();
    }
} // art_modern
} // labw