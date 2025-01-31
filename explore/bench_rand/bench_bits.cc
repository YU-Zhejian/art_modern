#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <mkl.h>

#include <gsl/gsl_rng.h>

#include <boost/random.hpp>

#include <absl/random/random.h>
#include <fstream>

#include "pcg-cpp-0.98/include/pcg_random.hpp"

#include "rprobs.hh"

class CustomRandomDevice {
public:
    using result_type = unsigned int;

    CustomRandomDevice()
        : rand_stream("/dev/random", std::ios::binary)
    {
        if (!rand_stream.is_open()) {
            throw std::system_error(errno, std::generic_category(), "Failed to open /dev/random");
        }
    }
    ~CustomRandomDevice() { rand_stream.close(); }

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return UINT_MAX; }

    result_type operator()()
    {
        result_type randomNumber;
        rand_stream.read(reinterpret_cast<char*>(&randomNumber), sizeof(randomNumber));
        if (rand_stream.fail()) {
            throw std::system_error(errno, std::generic_category(), "Failed to read from /dev/random");
        }

        return randomNumber;
    }

private:
    std::ifstream rand_stream;
};

class GslRandWrapper {
public:
    using result_type = unsigned long;
    GslRandWrapper(const gsl_rng_type* t) { r = gsl_rng_alloc(t); }
    ~GslRandWrapper() { gsl_rng_free(r); }
    result_type operator()() { return gsl_rng_get(r); }
    result_type min() { return gsl_rng_min(r); }
    result_type max() { return gsl_rng_max(r); }
    std::string name() const { return gsl_rng_name(r); }

private:
    gsl_rng* r;
};

namespace {
std::string formatWithCommas(std::size_t number)
{
    std::string numStr = std::to_string(number);
    int insertPosition = static_cast<int>(numStr.length()) - 3;

    while (insertPosition > 0) {
        numStr.insert(insertPosition, ",");
        insertPosition -= 3;
    }

    return numStr;
}

template <typename T> void bench_bits_stl(T& rng, const std::string& name)
{
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    std::vector<std::result_of_t<T()>> gen_bits(N_BASES);

    start = std::chrono::system_clock::now();
    for (int i = 0; i < N_TIMES; i++) {
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
    for (int i = 0; i < N_TIMES; i++) {
        viRngUniformBits(VSL_RNG_METHOD_UNIFORM_STD, stream, N_BASES, gen_bits.data());
    }
    end = std::chrono::system_clock::now();

    vslDeleteStream(&stream);
    std::cout << name << " (" << brng.NBits << " bits): "
              << formatWithCommas(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) << " us"
              << std::endl;
}

void bench_gsl(const gsl_rng_type* t)
{
    GslRandWrapper gsl_rand_wrapper { t };
    bench_bits_stl<GslRandWrapper>(gsl_rand_wrapper, "GSL::" + gsl_rand_wrapper.name());
}

void stl_main()
{
    CustomRandomDevice rng_custom_random_device;
    bench_bits_stl<CustomRandomDevice>(rng_custom_random_device, "CustomRandomDevice");

    std::mt19937 rng_mt19937 { seed() };
    bench_bits_stl<std::mt19937>(rng_mt19937, "std::mt19937");

    std::mt19937_64 rng_mt19937_64 { seed() };
    bench_bits_stl<std::mt19937_64>(rng_mt19937_64, "std::mt19937_64");

    std::ranlux48 rng_ranlux48 { seed() };
    bench_bits_stl<std::ranlux48>(rng_ranlux48, "std::ranlux48");

    std::ranlux24 rng_ranlux24 { seed() };
    bench_bits_stl<std::ranlux24>(rng_ranlux24, "std::ranlux24");

    std::ranlux48_base rng_ranlux48_base { seed() };
    bench_bits_stl<std::ranlux48_base>(rng_ranlux48_base, "std::ranlux48_base");

    std::ranlux24_base rng_ranlux24_base { seed() };
    bench_bits_stl<std::ranlux24_base>(rng_ranlux24_base, "std::ranlux24_base");

    std::knuth_b rng_knuth_b { seed() };
    bench_bits_stl<std::knuth_b>(rng_knuth_b, "std::knuth_b");

    std::minstd_rand0 rng_minstd_rand0 { seed() };
    bench_bits_stl<std::minstd_rand0>(rng_minstd_rand0, "std::minstd_rand0");

    std::minstd_rand rng_minstd_rand { seed() };
    bench_bits_stl<std::minstd_rand>(rng_minstd_rand, "std::minstd_rand");
    /**
     Very very slow
      std::random_device rng_random_device;
      std::vector<std::random_device::result_type> gen_bits_random_device(N_BASES);
      bench_bits_stl<std::random_device, std::random_device::result_type>(
          rng_random_device, gen_bits_random_device, "std::random_device");
    **/
}

void boost_main()
{
    boost::random::mt19937 rng_mt19937 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt19937>(rng_mt19937, "boost::random::mt19937");

    boost::random::mt19937_64 rng_mt19937_64 { seed() };
    bench_bits_stl<boost::random::mt19937_64>(rng_mt19937_64, "boost::random::mt19937_64");

    boost::random::ranlux48 rng_ranlux48 { seed() };
    bench_bits_stl<boost::random::ranlux48>(rng_ranlux48, "boost::random::ranlux48");

    boost::random::ranlux24 rng_ranlux24 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux24>(rng_ranlux24, "boost::random::ranlux24");

    boost::random::knuth_b rng_knuth_b { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::knuth_b>(rng_knuth_b, "boost::random::knuth_b");

    boost::random::minstd_rand0 rng_minstd_rand0 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::minstd_rand0>(rng_minstd_rand0, "boost::random::minstd_rand0");

    boost::random::minstd_rand rng_minstd_rand { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::minstd_rand>(rng_minstd_rand, "boost::random::minstd_rand");

    boost::random::ranlux48_base rng_ranlux48_base { seed() };
    bench_bits_stl<boost::random::ranlux48_base>(rng_ranlux48_base, "boost::random::ranlux48_base");

    boost::random::ranlux24_base rng_ranlux24_base { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux24_base>(rng_ranlux24_base, "boost::random::ranlux24_base");

    boost::random::rand48 rng_rand48 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::rand48>(rng_rand48, "boost::random::rand48");

    boost::random::ecuyer1988 rng_ecuyer1988 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ecuyer1988>(rng_ecuyer1988, "boost::random::ecuyer1988");

    boost::random::kreutzer1986 rng_kreutzer1986 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::kreutzer1986>(rng_kreutzer1986, "boost::random::kreutzer1986");

    boost::random::taus88 rng_taus88 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::taus88>(rng_taus88, "boost::random::taus88");

    boost::random::hellekalek1995 rng_hellekalek1995 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::hellekalek1995>(rng_hellekalek1995, "boost::random::hellekalek1995");

    boost::random::mt11213b rng_mt11213b { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt11213b>(rng_mt11213b, "boost::random::mt11213b");

    boost::random::lagged_fibonacci607 rng_lagged_fibonacci607 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci607>(rng_lagged_fibonacci607, "boost::random::lagged_fibonacci607");

    boost::random::lagged_fibonacci19937 rng_lagged_fibonacci19937 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci19937>(
        rng_lagged_fibonacci19937, "boost::random::lagged_fibonacci19937");

    boost::random::lagged_fibonacci9689 rng_lagged_fibonacci9689 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci9689>(
        rng_lagged_fibonacci9689, "boost::random::lagged_fibonacci9689");

    boost::random::lagged_fibonacci23209 rng_lagged_fibonacci23209 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci23209>(
        rng_lagged_fibonacci23209, "boost::random::lagged_fibonacci23209");

    boost::random::lagged_fibonacci1279 rng_lagged_fibonacci1279 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci1279>(
        rng_lagged_fibonacci1279, "boost::random::lagged_fibonacci1279");

    boost::random::lagged_fibonacci3217 rng_lagged_fibonacci3217 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci3217>(
        rng_lagged_fibonacci3217, "boost::random::lagged_fibonacci3217");

    boost::random::lagged_fibonacci4423 rng_lagged_fibonacci4423 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci4423>(
        rng_lagged_fibonacci4423, "boost::random::lagged_fibonacci4423");

    boost::random::lagged_fibonacci44497 rng_lagged_fibonacci44497 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::lagged_fibonacci44497>(
        rng_lagged_fibonacci44497, "boost::random::lagged_fibonacci44497");

    boost::random::ranlux3 rng_ranlux3 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux3>(rng_ranlux3, "boost::random::ranlux3");

    boost::random::ranlux4 rng_ranlux4 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux4>(rng_ranlux4, "boost::random::ranlux4");

    boost::random::ranlux64_3 rng_ranlux64_3 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_3>(rng_ranlux64_3, "boost::random::ranlux64_3");

    boost::random::ranlux64_4 rng_ranlux64_4 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_4>(rng_ranlux64_4, "boost::random::ranlux64_4");

    boost::random::ranlux3_01 rng_ranlux3_01 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux3_01>(rng_ranlux3_01, "boost::random::ranlux3_01");

    boost::random::ranlux4_01 rng_ranlux4_01 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux4_01>(rng_ranlux4_01, "boost::random::ranlux4_01");

    boost::random::ranlux64_3_01 rng_ranlux64_3_01 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_3_01>(rng_ranlux64_3_01, "boost::random::ranlux64_3_01");

    boost::random::ranlux64_4_01 rng_ranlux64_4_01 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_4_01>(rng_ranlux64_4_01, "boost::random::ranlux64_4_01");
}

void mkl_main()
{
    bench_bits_mkl(VSL_BRNG_MT19937, "MKL::VSL_BRNG_MT19937");
    bench_bits_mkl(VSL_BRNG_MT2203, "MKL::VSL_BRNG_MT2203");
    // This fails. No idea why.
    // bench_bits_mkl(VSL_BRNG_WH, "MKL::VSL_BRNG_WH");
    bench_bits_mkl(VSL_BRNG_SOBOL, "MKL::VSL_BRNG_SOBOL");
    bench_bits_mkl(VSL_BRNG_MCG31, "MKL::VSL_BRNG_MCG31");
    bench_bits_mkl(VSL_BRNG_R250, "MKL::VSL_BRNG_R250");
    bench_bits_mkl(VSL_BRNG_MRG32K3A, "MKL::VSL_BRNG_MRG32K3A");
    // This fails. No idea why.
    // bench_bits_mkl(VSL_BRNG_MCG59, "MKL::VSL_BRNG_MCG59");
    bench_bits_mkl(VSL_BRNG_NIEDERR, "MKL::VSL_BRNG_NIEDERR");
    bench_bits_mkl(VSL_BRNG_SFMT19937, "MKL::VSL_BRNG_SFMT19937");
    bench_bits_mkl(VSL_BRNG_ARS5, "MKL::VSL_BRNG_ARS5");
    bench_bits_mkl(VSL_BRNG_PHILOX4X32X10, "MKL::VSL_BRNG_PHILOX4X32X10");
    // These failed. No idea why.
    // bench_bits_mkl(VSL_BRNG_IABSTRACT, "MKL::VSL_BRNG_IABSTRACT");
    // bench_bits_mkl(VSL_BRNG_DABSTRACT, "MKL::VSL_BRNG_DABSTRACT");
    // bench_bits_mkl(VSL_BRNG_SABSTRACT, "MKL::VSL_BRNG_SABSTRACT");
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
    bench_gsl(gsl_rng_ranlxs0);
    bench_gsl(gsl_rng_ranlxs1);
    bench_gsl(gsl_rng_ranlxs2);
    bench_gsl(gsl_rng_ranlxd1);
    bench_gsl(gsl_rng_ranlxd2);
    bench_gsl(gsl_rng_ranlux);
    bench_gsl(gsl_rng_ranlux389);
    bench_gsl(gsl_rng_cmrg);
    bench_gsl(gsl_rng_mrg);
    bench_gsl(gsl_rng_taus);
    bench_gsl(gsl_rng_taus2);
    bench_gsl(gsl_rng_gfsr4);
    bench_gsl(gsl_rng_rand);
    bench_gsl(gsl_rng_random_bsd);
    bench_gsl(gsl_rng_random_libc5);
    bench_gsl(gsl_rng_random_glibc2);
    bench_gsl(gsl_rng_rand48);
    bench_gsl(gsl_rng_ranf);
    bench_gsl(gsl_rng_ranmar);
    bench_gsl(gsl_rng_r250);
    bench_gsl(gsl_rng_tt800);
    bench_gsl(gsl_rng_vax);
    bench_gsl(gsl_rng_transputer);
    bench_gsl(gsl_rng_randu);
    bench_gsl(gsl_rng_minstd);
    bench_gsl(gsl_rng_uni);
    bench_gsl(gsl_rng_uni32);
    bench_gsl(gsl_rng_slatec);
    bench_gsl(gsl_rng_zuf);
    bench_gsl(gsl_rng_knuthran2);
    bench_gsl(gsl_rng_knuthran2002);
    bench_gsl(gsl_rng_knuthran);
    bench_gsl(gsl_rng_borosh13);
    bench_gsl(gsl_rng_fishman18);
    bench_gsl(gsl_rng_fishman20);
    bench_gsl(gsl_rng_lecuyer21);
    bench_gsl(gsl_rng_waterman14);
    bench_gsl(gsl_rng_fishman2x);
    bench_gsl(gsl_rng_coveyou);
}
} // namespace

int main()
{
    stl_main();
    boost_main();
    mkl_main();
    absl_main();
    gsl_main();
    pcg_main();
    return EXIT_SUCCESS;
}
