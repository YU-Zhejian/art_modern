#include "gsl_rng_wrapper.hh"
#include "rprobs.hh"
#include "vigna.h"
#include "vmt19937_wrapper.hh"

// TestU01 library does not have C++ headers
extern "C" {
#include <TestU01.h>
}

#include <mkl.h>

#include <gsl/gsl_rng.h>

#include <boost/random.hpp>

#include <absl/random/random.h>

#include <pcg_random.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

namespace {
} // namespace

int main() noexcept
{
    std::string name = "CustomRandomDevice";

    char* rng_name = new char[name.size() + 1];
    std::strcpy(rng_name, name.c_str());
    unif01_Gen* gen = ufile_CreateReadBin("/dev/random", 4096);

    // Before running TestU01 battery
    auto* file = std::freopen("testu01_results.txt", "w", stdout);
    bbattery_SmallCrush(gen);
    std::fclose(file);
    unif01_DeleteExternGenBits(gen);
    return EXIT_SUCCESS;
}
