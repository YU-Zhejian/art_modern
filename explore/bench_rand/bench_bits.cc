/**
 *
 * @brief Benchmark various RNGs for generating random bits.
 *
 * Note: Those who failed TestU01 SmallCrush are NOT included.
 */
#include "gsl_rng_wrapper.hh"
#include "rprobs.hh"
#include "vigna.h"
#include "vmt19937_wrapper.hh"
#include "xoroshiro_wrapper.hh"

#include "arc4.hpp"
#include "gjrand.hpp"
#include "jsf.hpp"
#include "lehmer.hpp"
#include "sfc.hpp"
#include "splitmix.hpp"

#include <mkl.h>

#include <gsl/gsl_rng.h>

#include <boost/random.hpp>

#include <absl/base/attributes.h>
#include <absl/random/random.h>

#include <pcg_random.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t name_length = 65;

template <typename T> T geometric_mean(const std::vector<T>& data)
{
    double log_sum = 0.0;
    for (const auto& value : data) {
        log_sum += std::log(value);
    }
    return static_cast<T>(std::exp(log_sum / static_cast<double>(data.size())));
}

template <typename T> T mean(const std::vector<T>& data)
{
    return std::accumulate(data.begin(), data.end(), static_cast<T>(0)) / static_cast<T>(data.size());
}

template <typename T> T sd(const std::vector<T>& data, T mean_)
{
    T sum_squared_diff = 0;
    for (const auto& value : data) {
        T diff = value - mean_;
        sum_squared_diff += diff * diff;
    }
    return std::sqrt(sum_squared_diff / static_cast<T>(data.size() - 1));
}

std::string describe(const std::vector<std::size_t>& times)
{
    auto const mean_ = mean(times);
    std::ostringstream oss;
    oss << "gmean: " << std::setw(10) << formatWithCommas(geometric_mean(times)) << "; mean/sd: " << std::setw(15)
        << formatWithCommas(mean_) + "/" + formatWithCommas(sd(times, mean_));
    return oss.str();
}

template <typename T> void bench_bits_stl(T& rng, const std::string& name)
{
    std::vector<std::size_t> times {};

    std::vector<std::invoke_result_t<T>> gen_bits(N_BASES);

    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < N_TIMES; i++) {
            std::generate_n(gen_bits.begin(), N_BASES, [&rng]() { return rng(); });
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }
    std::string range;
    if (static_cast<std::int64_t>(rng.max()) - static_cast<std::int64_t>(rng.min()) + 1ULL == (1ULL << 32)) {
        range = "(32 bits)";
    } else if (static_cast<std::int64_t>(rng.max()) - static_cast<std::int64_t>(rng.min()) + 1ULL == (1ULL << 48)) {
        range = "(48 bits)";
    } else if (rng.min() == 0 && rng.max() == std::numeric_limits<std::uint64_t>::max()) {
        range = "(64 bits)";
    } else {
        range = "(" + std::to_string(rng.min()) + ", " + std::to_string(rng.max()) + ")";
    }

    std::cout << std::setw(name_length) << name + range + ": " << describe(times) << " us" << std::endl;
}

void bench_bits_mkl(const MKL_INT type, const std::string& name)
{
    VSLStreamStatePtr stream = nullptr;
    vslNewStream(&stream, type, seed());
    VSLBRngProperties brng;
    vslGetBrngProperties(type, &brng);
    std::vector<std::size_t> times {};
    std::vector<std::uint32_t> gen_bits {};
    gen_bits.resize(N_BASES);

    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < N_TIMES; i++) {
            viRngUniformBits(VSL_RNG_METHOD_UNIFORM_STD, stream, N_BASES, gen_bits.data());
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }

    vslDeleteStream(&stream);
    std::cout << std::setw(name_length) << name + "(" + std::to_string(brng.NBits) + " bits): " << describe(times)
              << " us" << std::endl;
}

template <typename VMT19937BulkRandomDeviceImpl> void bench_bits_vmt19937(const std::string& name)
{
    VMT19937BulkRandomDeviceImpl rng {};
    std::vector<std::size_t> times {};
    std::vector<std::uint32_t> gen_bits {};
    gen_bits.resize(N_BASES);

    for (std::size_t j = 0; j < N_REPLICA; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < N_TIMES; i++) {
            rng.gen(gen_bits);
        }
        auto end = std::chrono::high_resolution_clock::now();
        times.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
    }

    std::string range;
    if (static_cast<std::int64_t>(rng.max()) - static_cast<std::int64_t>(rng.min()) + 1ULL == (1ULL << 32)) {
        range = "(32 bits)";
    } else if (static_cast<std::int64_t>(rng.max()) - static_cast<std::int64_t>(rng.min()) + 1ULL == (1ULL << 48)) {
        range = "(48 bits)";
    } else if (rng.min() == 0 && rng.max() == std::numeric_limits<std::uint64_t>::max()) {
        range = "(64 bits)";
    } else {
        range = "(" + std::to_string(rng.min()) + ", " + std::to_string(rng.max()) + ")";
    }

    std::cout << std::setw(name_length) << name + range + ": " << describe(times) << " us" << std::endl;
}

ABSL_ATTRIBUTE_ALWAYS_INLINE void bench_gsl(const gsl_rng_type* t)
{
    GslRngWrapper gsl_rand_wrapper { t };
    bench_bits_stl<GslRngWrapper>(gsl_rand_wrapper, "GSL::" + gsl_rand_wrapper.name());
}

[[maybe_unused]] void stl_main()
{
    CustomRandomDevice rng_custom_random_device;
    bench_bits_stl<CustomRandomDevice>(rng_custom_random_device, "CustomRandomDevice");

    DumbRandomDevice rng_dumb_random_device;
    bench_bits_stl<DumbRandomDevice>(rng_dumb_random_device, "DumbRandomDevice");

    std::mt19937 rng_mt19937 { seed() };
    bench_bits_stl<std::mt19937>(rng_mt19937, "std::mt19937");

    std::mt19937_64 rng_mt19937_64 { seed() };
    bench_bits_stl<std::mt19937_64>(rng_mt19937_64, "std::mt19937_64");

    std::ranlux48 rng_ranlux48 { seed() };
    bench_bits_stl<std::ranlux48>(rng_ranlux48, "std::ranlux48");
}

[[maybe_unused]] void boost_main()
{
    boost::random::mt19937 rng_mt19937 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt19937>(rng_mt19937, "boost::random::mt19937");

    boost::random::mt19937_64 rng_mt19937_64 { seed() };
    bench_bits_stl<boost::random::mt19937_64>(rng_mt19937_64, "boost::random::mt19937_64");

    boost::random::ranlux48 rng_ranlux48 { seed() };
    bench_bits_stl<boost::random::ranlux48>(rng_ranlux48, "boost::random::ranlux48");

    boost::random::taus88 rng_taus88 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::taus88>(rng_taus88, "boost::random::taus88");

    boost::random::mt11213b rng_mt11213b { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt11213b>(rng_mt11213b, "boost::random::mt11213b");

    boost::random::ranlux64_3 rng_ranlux64_3 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_3>(rng_ranlux64_3, "boost::random::ranlux64_3");

    boost::random::ranlux64_4 rng_ranlux64_4 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_4>(rng_ranlux64_4, "boost::random::ranlux64_4");
}

[[maybe_unused]] void mkl_main()
{
    bench_bits_mkl(VSL_BRNG_MT19937, "MKL::VSL_BRNG_MT19937");
    bench_bits_mkl(VSL_BRNG_MT2203, "MKL::VSL_BRNG_MT2203");
    bench_bits_mkl(VSL_BRNG_SFMT19937, "MKL::VSL_BRNG_SFMT19937");
    bench_bits_mkl(VSL_BRNG_ARS5, "MKL::VSL_BRNG_ARS5");
    bench_bits_mkl(VSL_BRNG_PHILOX4X32X10, "MKL::VSL_BRNG_PHILOX4X32X10");
    // Yet another /dev/random
    // bench_bits_mkl(VSL_BRNG_NONDETERM, "MKL::VSL_BRNG_NONDETERM");
}

[[maybe_unused]] void absl_main()
{
    absl::BitGen rng_bitgen {};
    bench_bits_stl<absl::BitGen>(rng_bitgen, "absl::BitGen");

    absl::InsecureBitGen rng_insecure_bitgen {};
    bench_bits_stl<absl::InsecureBitGen>(rng_insecure_bitgen, "absl::InsecureBitGen");
}

[[maybe_unused]] void pcg_main()
{
    pcg32 rng_pcg32 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg32>(rng_pcg32, "PCG::pcg32");

    pcg64 rng_pcg64 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg64>(rng_pcg64, "PCG::pcg64");

    pcg32_fast rng_pcg32_fast { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg32_fast>(rng_pcg32_fast, "PCG::pcg32_fast");

    pcg64_fast rng_pcg64_fast { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg64_fast>(rng_pcg64_fast, "PCG::pcg64_fast");

    pcg32_oneseq_once_insecure rng_pcg32_oneseq_once_insecure { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg32_oneseq_once_insecure>(rng_pcg32_oneseq_once_insecure, "PCG::pcg32_oneseq_once_insecure");

    pcg64_oneseq_once_insecure rng_pcg64_oneseq_once_insecure { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg64_oneseq_once_insecure>(rng_pcg64_oneseq_once_insecure, "PCG::pcg64_oneseq_once_insecure");
}

[[maybe_unused]] void gsl_main()
{
    bench_gsl(gsl_rng_mt19937);
    bench_gsl(gsl_rng_mt19937_1999);
    bench_gsl(gsl_rng_mt19937_1998);
    bench_gsl(gsl_rng_ranlxd1);
    bench_gsl(gsl_rng_ranlxd2);
    bench_gsl(gsl_rng_taus);
    bench_gsl(gsl_rng_taus2);
    bench_gsl(gsl_rng_gfsr4);
}

[[maybe_unused]] void vmt19937_main()
{
    VMT19937RandomDevice rng_vmt19937_random_device {};
    bench_bits_stl<decltype(rng_vmt19937_random_device)>(rng_vmt19937_random_device, "VMT19937RandomDevice");

    VSFMT19937RandomDevice rng_vsfmt19937_random_device {};
    bench_bits_stl<decltype(rng_vsfmt19937_random_device)>(rng_vsfmt19937_random_device, "VSFMT19937RandomDevice");

    bench_bits_vmt19937<VMT19937BulkRandomDevice>("VMT19937BulkRandomDevice");
    bench_bits_vmt19937<VSFMT19937BulkRandomDevice>("VSFMT19937BulkRandomDevice");
}

[[maybe_unused]] void xso_main()
{
    XoroshiroWrapper<old::xoroshiro_2x32_star, uint32_t> x01 {};
    bench_bits_stl<decltype(x01)>(x01, "xoroshiro::2x32*");

    XoroshiroWrapper<old::xoshiro_4x32_plus, uint32_t, 4> x03 {};
    bench_bits_stl<decltype(x03)>(x03, "xoshiro::4x32+");

    XoroshiroWrapper<old::xoshiro_4x32_plus_plus, uint32_t, 4> x04 {};
    bench_bits_stl<decltype(x04)>(x04, "xoshiro::4x32++");

    XoroshiroWrapper<old::xoshiro_4x32_star_star, uint32_t, 4> x05 {};
    bench_bits_stl<decltype(x05)>(x05, "xoshiro::4x32**");

    XoroshiroWrapper<old::xoroshiro_2x64_plus, uint64_t, 2> x06 {};
    bench_bits_stl<decltype(x06)>(x06, "xoroshiro::2x64+");

    XoroshiroWrapper<old::xoroshiro_2x64_plus_plus, uint64_t, 2> x07 {};
    bench_bits_stl<decltype(x07)>(x07, "xoroshiro::2x64++");

    XoroshiroWrapper<old::xoroshiro_2x64_star_star, uint64_t, 2> x08 {};
    bench_bits_stl<decltype(x08)>(x08, "xoroshiro::2x64**");

    XoroshiroWrapper<old::xoshiro_4x64_plus, uint64_t, 4> x09 {};
    bench_bits_stl<decltype(x09)>(x09, "xoshiro::4x64+");

    XoroshiroWrapper<old::xoshiro_4x64_plus_plus, uint64_t, 4> x10 {};
    bench_bits_stl<decltype(x10)>(x10, "xoshiro::4x64++");

    XoroshiroWrapper<old::xoshiro_4x64_star_star, uint64_t, 4> x11 {};
    bench_bits_stl<decltype(x11)>(x11, "xoshiro::4x64**");

    XoroshiroWrapper<old::xoshiro_8x64_plus, uint64_t, 8> x12 {};
    bench_bits_stl<decltype(x12)>(x12, "xoshiro::8x64+");

    XoroshiroWrapper<old::xoshiro_8x64_plus_plus, uint64_t, 8> x13 {};
    bench_bits_stl<decltype(x13)>(x13, "xoshiro::8x64++");

    XoroshiroWrapper<old::xoshiro_8x64_star_star, uint64_t, 8> x14 {};
    bench_bits_stl<decltype(x14)>(x14, "xoshiro::8x64**");

    XoroshiroWrapper<old::xoroshiro_16x64_star, uint64_t, 16> x15 {};
    bench_bits_stl<decltype(x15)>(x15, "xoroshiro::16x64*");

    XoroshiroWrapper<old::xoroshiro_16x64_star_star, uint64_t, 16> x16 {};
    bench_bits_stl<decltype(x16)>(x16, "xoroshiro::16x64**");

    XoroshiroWrapper<old::xoroshiro_16x64_plus_plus, uint64_t, 16> x17 {};
    bench_bits_stl<decltype(x17)>(x17, "xoroshiro::16x64++");
}

[[maybe_unused]] void other_rngs_main()
{
    arc4_rand32 rng_arc4 { static_cast<arc4_rand32::result_type>(seed()) };
    bench_bits_stl<arc4_rand32>(rng_arc4, "others::arc4_rand32");

    arc4_rand64 rng_arc4_64 { static_cast<arc4_rand64::result_type>(seed()) };
    bench_bits_stl<arc4_rand64>(rng_arc4_64, "others::arc4_rand64");

    gjrand32 rng_gjrand32 { static_cast<gjrand32::result_type>(seed()) };
    bench_bits_stl<gjrand32>(rng_gjrand32, "others::gjrand32");

    gjrand64 rng_gjrand63 { static_cast<gjrand64::result_type>(seed()) };
    bench_bits_stl<gjrand64>(rng_gjrand63, "others::gjrand64");

    jsf32 rng_jsf32 { static_cast<jsf32::result_type>(seed()) };
    bench_bits_stl<jsf32>(rng_jsf32, "others::jsf32");

    mcg128 rng_mcg128 { static_cast<mcg128::result_type>(seed()) };
    bench_bits_stl<mcg128>(rng_mcg128, "others::mcg128");

    mcg128_fast rng_mcg128_fast { static_cast<mcg128_fast::result_type>(seed()) };
    bench_bits_stl<mcg128_fast>(rng_mcg128_fast, "others::mcg128_fast");

    sfc32 rng_sfc32 { static_cast<sfc32::result_type>(seed()) };
    bench_bits_stl<sfc32>(rng_sfc32, "others::sfc32");

    sfc64 rng_sfc64 { static_cast<sfc64::result_type>(seed()) };
    bench_bits_stl<sfc64>(rng_sfc64, "others::sfc64");

    splitmix32 rng_splitmix32 { static_cast<splitmix32::result_type>(seed()) };
    bench_bits_stl<splitmix32>(rng_splitmix32, "others::splitmix32");

    splitmix64 rng_splitmix64 { static_cast<splitmix64::result_type>(seed()) };
    bench_bits_stl<splitmix64>(rng_splitmix64, "others::splitmix64");
}

} // namespace

int main() noexcept
{
    stl_main();
    boost_main();
    mkl_main();
    absl_main();
    gsl_main();
    pcg_main();
    xso_main();
    vmt19937_main();
    other_rngs_main();
    return EXIT_SUCCESS;
}
