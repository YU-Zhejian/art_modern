/**
 *
 * Note: Those who failed TestU01 are NOT included.
 */
#include "gsl_rng_wrapper.hh"
#include "rprobs.hh"
#include "vigna.h"
#include "vmt19937_wrapper.hh"
#include "xoroshiro_wrapper.hh"

#include <mkl.h>

#include <gsl/gsl_rng.h>

#include <boost/random.hpp>

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
template <typename T> void bench_bits_stl(T& rng, const std::string& name)
{
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    std::vector<std::invoke_result_t<T>> gen_bits(N_BASES);

    start = std::chrono::system_clock::now();
    for (std::size_t i = 0; i < N_TIMES; i++) {
        std::generate_n(gen_bits.begin(), N_BASES, [&rng]() { return rng(); });
    }
    end = std::chrono::system_clock::now();

    std::cout << name << "(" << std::to_string(rng.min()) << ", " << std::to_string(rng.max())
              << "): " << formatWithCommas(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())
              << " us" << std::endl;
}

void bench_bits_mkl(const MKL_INT type, const std::string& name)
{
    VSLStreamStatePtr stream = nullptr;
    vslNewStream(&stream, type, seed());
    VSLBRngProperties brng;
    vslGetBrngProperties(type, &brng);
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    std::vector<std::uint32_t> gen_bits {};
    gen_bits.resize(N_BASES);

    start = std::chrono::system_clock::now();
    for (std::size_t i = 0; i < N_TIMES; i++) {
        viRngUniformBits(VSL_RNG_METHOD_UNIFORM_STD, stream, N_BASES, gen_bits.data());
    }
    end = std::chrono::system_clock::now();

    vslDeleteStream(&stream);
    std::cout << name << " (" << brng.NBits << " bits): "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) << " us"
              << std::endl;
}

template <typename VMT19937BulkRandomDeviceImpl> void bench_bits_vmt19937(const std::string& name)
{
    VMT19937BulkRandomDeviceImpl rng {};
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    std::vector<std::uint32_t> gen_bits {};
    gen_bits.resize(N_BASES);

    start = std::chrono::system_clock::now();
    for (std::size_t i = 0; i < N_TIMES; i++) {
        rng.gen(gen_bits);
    }
    end = std::chrono::system_clock::now();

    std::cout << name << "(" << std::to_string(rng.min()) << ", " << std::to_string(rng.max())
              << "): " << formatWithCommas(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())
              << " us" << std::endl;
}

void bench_gsl(const gsl_rng_type* t)
{
    GslRngWrapper gsl_rand_wrapper { t };
    bench_bits_stl<GslRngWrapper>(gsl_rand_wrapper, "GSL::" + gsl_rand_wrapper.name());
}

void stl_main()
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

    std::ranlux24 rng_ranlux24 { seed() };
    bench_bits_stl<std::ranlux24>(rng_ranlux24, "std::ranlux24");

    std::knuth_b rng_knuth_b { seed() };
    bench_bits_stl<std::knuth_b>(rng_knuth_b, "std::knuth_b");
}

void boost_main()
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

void mkl_main()
{
    bench_bits_mkl(VSL_BRNG_MT19937, "MKL::VSL_BRNG_MT19937");
    bench_bits_mkl(VSL_BRNG_MT2203, "MKL::VSL_BRNG_MT2203");
    bench_bits_mkl(VSL_BRNG_SFMT19937, "MKL::VSL_BRNG_SFMT19937");
    bench_bits_mkl(VSL_BRNG_ARS5, "MKL::VSL_BRNG_ARS5");
    bench_bits_mkl(VSL_BRNG_PHILOX4X32X10, "MKL::VSL_BRNG_PHILOX4X32X10");
    bench_bits_mkl(VSL_BRNG_NONDETERM, "MKL::VSL_BRNG_NONDETERM");
}

void absl_main()
{
    absl::BitGen rng_bitgen {};
    bench_bits_stl<absl::BitGen>(rng_bitgen, "absl::BitGen");

    absl::InsecureBitGen rng_insecure_bitgen {};
    bench_bits_stl<absl::InsecureBitGen>(rng_insecure_bitgen, "absl::InsecureBitGen");
}

void pcg_main()
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

void gsl_main()
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

void mt19937_main()
{
    std::mt19937 rng_stl_mt19937 { seed() };
    bench_bits_stl<std::mt19937>(rng_stl_mt19937, "std::mt19937");

    boost::random::mt19937 rng_boost_mt19937 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt19937>(rng_boost_mt19937, "boost::random::mt19937");

    bench_bits_mkl(VSL_BRNG_MT19937, "MKL::VSL_BRNG_MT19937");
    bench_bits_mkl(VSL_BRNG_SFMT19937, "MKL::VSL_BRNG_SFMT19937");
    bench_gsl(gsl_rng_mt19937);

    bench_bits_vmt19937<VMT19937BulkRandomDevice>("VMT19937BulkRandomDevice");
    bench_bits_vmt19937<VSFMT19937BulkRandomDevice>("VSFMT19937BulkRandomDevice");
}

void xoshiro_main()
{
    XoroshiroWrapper<old::xoroshiro_2x32_star, uint32_t> x01 {};
    bench_bits_stl<decltype(x01)>(x01, "xoroshiro::2x32*");

    XoroshiroWrapper<old::xoroshiro_2x32_star_star, uint32_t> x02 {};
    bench_bits_stl<decltype(x02)>(x02, "xoroshiro::2x32**");

    XoroshiroWrapper<old::xoshiro_4x32_plus, uint32_t, 4> x03 {};
    bench_bits_stl<decltype(x03)>(x03, "xoshiro::4x32+");

    XoroshiroWrapper<old::xoshiro_4x32_star_star, uint32_t, 4> x04 {};
    bench_bits_stl<decltype(x04)>(x04, "xoshiro::4x32**");

    XoroshiroWrapper<old::xoshiro_4x32_plus_plus, uint32_t, 4> x05 {};
    bench_bits_stl<decltype(x05)>(x05, "xoshiro::4x32++");
}

} // namespace

int main() noexcept
{
    mt19937_main();
    stl_main();
    boost_main();
    mkl_main();
    absl_main();
    gsl_main();
    pcg_main();
    xoshiro_main();
    return EXIT_SUCCESS;
}
