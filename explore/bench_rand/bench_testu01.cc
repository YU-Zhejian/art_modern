#include "class_utils.hh"
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

#include <boost/filesystem.hpp>
#include <boost/process.hpp>
#include <boost/random.hpp>

#include <absl/random/random.h>

#include <pcg_random.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <utility>
#include <vector>

static std::string testu01_path;

namespace bp = boost::process;

class ContinuousDataWriter {

public:
    ContinuousDataWriter(bp::async_pipe& pipe, bp::child& proc)
        : stdin_pipe_(pipe)
        , process_(proc)
    {
    }
    ~ContinuousDataWriter() { stop(); }
    DELETE_COPY_MOVE(ContinuousDataWriter)

    void write(const char* data_to_send = nullptr, std::size_t size = 0)
    {
        if (good()) {
            try {
                actual_bytes_transferred_ += boost::asio::write(stdin_pipe_, boost::asio::buffer(data_to_send, size));
            } catch (const std::exception& e) {
                // std::cerr << "Write error: " << e.what() << std::endl;
                stop();
            }
        }
    }
    void stop() { should_continue_ = false; }
    [[nodiscard]] std::size_t actual_bytes_transferred() const { return actual_bytes_transferred_; }
    [[nodiscard]] bool good() const { return should_continue_ && process_.running() && stdin_pipe_.is_open(); }

private:
    bp::async_pipe& stdin_pipe_;
    bool should_continue_ = true;
    bp::child& process_;
    std::size_t actual_bytes_transferred_ = 0;
};

class BenchmarkingHelper {
private:
    std::string name_;

public:
    explicit BenchmarkingHelper(std::string name_)
        : name_(std::move(name_))
    {
    }
    virtual ~BenchmarkingHelper() = default;
    DELETE_COPY_MOVE(BenchmarkingHelper)
    void run()
    {
        std::ofstream f { "testu01_results.txt", std::ios::app };

        f << ">" << name_ << "\n";

        boost::asio::io_context io_context;
        bp::async_pipe stdin_pipe { io_context };
        std::future<std::string> stdout_buffer;
        // Create the child process
        bp::child c(
            testu01_path, bp::std_in<stdin_pipe, bp::std_out> stdout_buffer, bp::std_err > bp::null, io_context);
        const auto range = this->range();
        std::cerr << name_ << ": Ranged " << formatWithCommas(range) << "." << std::endl;
        if (range < std::numeric_limits<std::uint32_t>::max()) {
            std::cerr << name_ << ": Range too small, skipping." << std::endl;
            f << "Range too small, skipping.\n";
            f.close();
            return;
        }
        std::cerr << name_ << ": Subprocess started." << std::endl;

        ContinuousDataWriter cdw { stdin_pipe, c };
        std::thread io_thread([&io_context]() { io_context.run(); });
        auto rand_bytes_adaptor_thread = create_rand_bytes_adaptor_thread(cdw);
        rand_bytes_adaptor_thread.join();
        stdin_pipe.close();
        c.wait();
        io_thread.join();
        auto stdout_contents = stdout_buffer.get();
        f.write(stdout_contents.c_str(), stdout_contents.size());
        f.close();

        std::cerr << name_ << ": Subprocess finished with " << formatWithCommas(cdw.actual_bytes_transferred())
                  << " bytes." << std::endl;
    }

    virtual std::thread create_rand_bytes_adaptor_thread(ContinuousDataWriter& cdw) = 0;
    virtual std::size_t range() = 0;
};

template <typename T, typename N = std::uint32_t> class StlBenchmarkingHelper : public BenchmarkingHelper {
private:
    T& rng_;

public:
    StlBenchmarkingHelper(T& rng, std::string name_)
        : BenchmarkingHelper(std::move(name_))
        , rng_(rng)
    {
    }
    ~StlBenchmarkingHelper() override = default;
    DELETE_COPY_MOVE(StlBenchmarkingHelper)
    std::thread create_rand_bytes_adaptor_thread(ContinuousDataWriter& cdw) override
    {
        return std::thread([&cdw, this]() { rand_bytes_adaptor(cdw); });
    }

    void rand_bytes_adaptor(ContinuousDataWriter& cdw)
    {
        auto uid = std::uniform_int_distribution<N>(0, std::numeric_limits<N>::max());
        std::vector<N> buffer { };
        buffer.resize(4096);
        while (cdw.good()) {
            std::generate_n(buffer.begin(), buffer.size(), [&uid, this]() { return uid(rng_); });
            cdw.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(N));
        }
    }
    std::size_t range() override { return static_cast<std::size_t>(rng_.max() - rng_.min()); }
};

template <typename N = std::uint32_t> class MklBenchmarkingHelper : public BenchmarkingHelper {
private:
    VSLStreamStatePtr stream_ = nullptr;
    MKL_INT type_;

public:
    MklBenchmarkingHelper(const MKL_INT type, std::string name_)
        : BenchmarkingHelper(std::move(name_))
        , type_(type)
    {
        vslNewStream(&stream_, type_, seed());
    }
    ~MklBenchmarkingHelper() override { vslDeleteStream(&stream_); };
    DELETE_COPY_MOVE(MklBenchmarkingHelper)
    std::thread create_rand_bytes_adaptor_thread(ContinuousDataWriter& cdw) override
    {
        return std::thread([&cdw, this]() { rand_bytes_adaptor(cdw); });
    }

    void rand_bytes_adaptor(ContinuousDataWriter& cdw)
    {
        auto uid = std::uniform_int_distribution<N>(0, std::numeric_limits<N>::max());
        std::vector<N> buffer { };
        buffer.resize(4096);
        while (cdw.good()) {
            viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream_, buffer.size(), buffer.data());
            cdw.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(N));
        }
    }
    std::size_t range() override
    {
        VSLBRngProperties brng;
        vslGetBrngProperties(type_, &brng);
        return static_cast<std::size_t>(1ULL << brng.NBits);
    }
};

template <typename N = std::uint32_t> class GslBenchmarkingHelper : public BenchmarkingHelper {
private:
    gsl_rng* gen_;

public:
    GslBenchmarkingHelper(const gsl_rng_type* type, std::string name_)
        : BenchmarkingHelper(std::move(name_))
        , gen_(gsl_rng_alloc(type))
    {
        gsl_rng_set(gen_, seed());
    }
    ~GslBenchmarkingHelper() override { gsl_rng_free(gen_); }
    DELETE_COPY_MOVE(GslBenchmarkingHelper)
    std::thread create_rand_bytes_adaptor_thread(ContinuousDataWriter& cdw) override
    {
        return std::thread([&cdw, this]() { rand_bytes_adaptor(cdw); });
    }

    void rand_bytes_adaptor(ContinuousDataWriter& cdw)
    {
        std::vector<N> buffer { };
        buffer.resize(4096);
        while (cdw.good()) {
            std::generate_n(buffer.begin(), buffer.size(),
                [this]() { return gsl_rng_uniform_int(gen_, std::numeric_limits<N>::max()); });

            cdw.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(N));
        }
    }
    std::size_t range() override { return static_cast<std::size_t>(gsl_rng_max(gen_) - gsl_rng_min(gen_)); }
};

namespace {

template <typename T> void bench_bits_stl(T& rng, const std::string& name)
{
    StlBenchmarkingHelper { rng, name }.run();
}

void bench_gsl(const gsl_rng_type* t)
{
    auto* rng = gsl_rng_alloc(t);
    GslBenchmarkingHelper { t, "GSL::" + std::string(gsl_rng_name(rng)) }.run();
    gsl_rng_free(rng);
}

void bench_bits_mkl(const MKL_INT type, const std::string& name) { MklBenchmarkingHelper { type, name }.run(); }

[[maybe_unused]] void stl_main()
{
    CustomRandomDevice rng_custom_random_device;
    bench_bits_stl<decltype(rng_custom_random_device)>(rng_custom_random_device, "CustomRandomDevice");
    DumbRandomDevice rng_dumb_random_device;
    bench_bits_stl<decltype(rng_dumb_random_device)>(rng_dumb_random_device, "DumbRandomDevice");

    // All tests were passed
    std::mt19937 rng_mt19937 { seed() };
    bench_bits_stl<decltype(rng_mt19937)>(rng_mt19937, "std::mt19937");

    // All tests were passed
    std::mt19937_64 rng_mt19937_64 { seed() };
    bench_bits_stl<std::mt19937_64>(rng_mt19937_64, "std::mt19937_64");

    // All tests were passed
    std::ranlux48 rng_ranlux48 { seed() };
    bench_bits_stl<std::ranlux48>(rng_ranlux48, "std::ranlux48");

    // Range too small, skipping.
    // std::ranlux24 rng_ranlux24 { seed() };
    // bench_bits_stl<std::ranlux24>(rng_ranlux24, "std::ranlux24");

    //        Test                          p-value
    // ----------------------------------------------
    //  4  SimpPoker                        eps
    //  7  WeightDistrib                  1.1e-16
    // ----------------------------------------------
    // All other tests were passed
    // std::ranlux48_base rng_ranlux48_base { seed() };
    // bench_bits_stl<std::ranlux48_base>(rng_ranlux48_base, "std::ranlux48_base");

    //        Test                          p-value
    // ----------------------------------------------
    //  4  SimpPoker                        eps
    //  5  CouponCollector                 3.6e-4
    //  7  WeightDistrib                    eps
    // ----------------------------------------------
    // All other tests were passed
    // std::ranlux24_base rng_ranlux24_base { seed() };
    // bench_bits_stl<std::ranlux24_base>(rng_ranlux24_base, "std::ranlux24_base");

    // Range too small, skipping.
    // std::knuth_b rng_knuth_b { seed() };
    // bench_bits_stl<std::knuth_b>(rng_knuth_b, "std::knuth_b");

    //        Test                          p-value
    // ----------------------------------------------
    //  1  BirthdaySpacings                 eps
    //  2  Collision                      1 - eps1
    // ----------------------------------------------
    // All other tests were passed
    // std::minstd_rand0 rng_minstd_rand0 { seed() };
    // bench_bits_stl<std::minstd_rand0>(rng_minstd_rand0, "std::minstd_rand0");

    //        Test                          p-value
    // ----------------------------------------------
    //  1  BirthdaySpacings                 eps
    //  2  Collision                      1 - eps1
    // ----------------------------------------------
    // All other tests were passed
    // std::minstd_rand rng_minstd_rand { seed() };
    // bench_bits_stl<std::minstd_rand>(rng_minstd_rand, "std::minstd_rand");
}

[[maybe_unused]] void boost_main()
{
    // All tests were passed
    boost::random::mt19937 rng_mt19937 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt19937>(rng_mt19937, "boost::random::mt19937");

    // All tests were passed
    boost::random::mt19937_64 rng_mt19937_64 { seed() };
    bench_bits_stl<boost::random::mt19937_64>(rng_mt19937_64, "boost::random::mt19937_64");

    // All tests were passed
    boost::random::ranlux48 rng_ranlux48 { seed() };
    bench_bits_stl<boost::random::ranlux48>(rng_ranlux48, "boost::random::ranlux48");

    // Range too small, skipping.
    // boost::random::ranlux24 rng_ranlux24 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux24>(rng_ranlux24, "boost::random::ranlux24");

    // Range too small, skipping.
    // boost::random::knuth_b rng_knuth_b { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::knuth_b>(rng_knuth_b, "boost::random::knuth_b");

    // Range too small, skipping.
    // boost::random::minstd_rand0 rng_minstd_rand0 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::minstd_rand0>(rng_minstd_rand0, "boost::random::minstd_rand0");

    // Range too small, skipping.
    // boost::random::minstd_rand rng_minstd_rand { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::minstd_rand>(rng_minstd_rand, "boost::random::minstd_rand");

    //        Test                          p-value
    // ----------------------------------------------
    //  4  SimpPoker                        eps
    //  7  WeightDistrib                    eps
    // ----------------------------------------------
    // All other tests were passed
    // boost::random::ranlux48_base rng_ranlux48_base { seed() };
    // bench_bits_stl<boost::random::ranlux48_base>(rng_ranlux48_base, "boost::random::ranlux48_base");

    // Range too small, skipping.
    // boost::random::ranlux24_base rng_ranlux24_base { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux24_base>(rng_ranlux24_base, "boost::random::ranlux24_base");

    //        Test                          p-value
    // ----------------------------------------------
    //  1  BirthdaySpacings                 eps
    //  2  Collision                      1 - eps1
    // ----------------------------------------------
    // All other tests were passed
    // boost::random::rand48 rng_rand48 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::rand48>(rng_rand48, "boost::random::rand48");

    // Range too small, skipping.
    // boost::random::ecuyer1988 rng_ecuyer1988 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ecuyer1988>(rng_ecuyer1988, "boost::random::ecuyer1988");

    // Range too small, skipping.
    // boost::random::kreutzer1986 rng_kreutzer1986 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::kreutzer1986>(rng_kreutzer1986, "boost::random::kreutzer1986");

    // All tests were passed
    boost::random::taus88 rng_taus88 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::taus88>(rng_taus88, "boost::random::taus88");

    // Range too small, skipping.
    // boost::random::hellekalek1995 rng_hellekalek1995 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::hellekalek1995>(rng_hellekalek1995, "boost::random::hellekalek1995");

    // All tests were passed
    boost::random::mt11213b rng_mt11213b { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::mt11213b>(rng_mt11213b, "boost::random::mt11213b");

    // Range too small, skipping.
    // boost::random::lagged_fibonacci607 rng_lagged_fibonacci607 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci607>(rng_lagged_fibonacci607,
    // "boost::random::lagged_fibonacci607"); boost::random::lagged_fibonacci19937 rng_lagged_fibonacci19937 {
    // static_cast<unsigned int>(seed()) }; bench_bits_stl<boost::random::lagged_fibonacci19937>(
    //     rng_lagged_fibonacci19937, "boost::random::lagged_fibonacci19937");
    // boost::random::lagged_fibonacci9689 rng_lagged_fibonacci9689 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci9689>(
    //     rng_lagged_fibonacci9689, "boost::random::lagged_fibonacci9689");
    // boost::random::lagged_fibonacci23209 rng_lagged_fibonacci23209 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci23209>(
    //     rng_lagged_fibonacci23209, "boost::random::lagged_fibonacci23209");
    // boost::random::lagged_fibonacci1279 rng_lagged_fibonacci1279 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci1279>(
    //     rng_lagged_fibonacci1279, "boost::random::lagged_fibonacci1279");
    // boost::random::lagged_fibonacci3217 rng_lagged_fibonacci3217 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci3217>(
    //     rng_lagged_fibonacci3217, "boost::random::lagged_fibonacci3217");
    // boost::random::lagged_fibonacci4423 rng_lagged_fibonacci4423 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci4423>(
    //     rng_lagged_fibonacci4423, "boost::random::lagged_fibonacci4423");
    // boost::random::lagged_fibonacci44497 rng_lagged_fibonacci44497 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::lagged_fibonacci44497>(
    //     rng_lagged_fibonacci44497, "boost::random::lagged_fibonacci44497");

    // Range too small, skipping.
    // boost::random::ranlux3 rng_ranlux3 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux3>(rng_ranlux3, "boost::random::ranlux3");

    // Range too small, skipping.
    // boost::random::ranlux4 rng_ranlux4 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux4>(rng_ranlux4, "boost::random::ranlux4");

    // All tests were passed
    boost::random::ranlux64_3 rng_ranlux64_3 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_3>(rng_ranlux64_3, "boost::random::ranlux64_3");

    // All tests were passed
    boost::random::ranlux64_4 rng_ranlux64_4 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<boost::random::ranlux64_4>(rng_ranlux64_4, "boost::random::ranlux64_4");

    // Range too small, skipping.
    // boost::random::ranlux3_01 rng_ranlux3_01 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux3_01>(rng_ranlux3_01, "boost::random::ranlux3_01");
    // boost::random::ranlux4_01 rng_ranlux4_01 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux4_01>(rng_ranlux4_01, "boost::random::ranlux4_01");
    // boost::random::ranlux64_3_01 rng_ranlux64_3_01 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux64_3_01>(rng_ranlux64_3_01, "boost::random::ranlux64_3_01");
    // boost::random::ranlux64_4_01 rng_ranlux64_4_01 { static_cast<unsigned int>(seed()) };
    // bench_bits_stl<boost::random::ranlux64_4_01>(rng_ranlux64_4_01, "boost::random::ranlux64_4_01");
}

[[maybe_unused]] void absl_main()
{
    // All tests were passed
    absl::BitGen rng_bitgen { };
    bench_bits_stl<absl::BitGen>(rng_bitgen, "absl::BitGen");

    // All tests were passed
    absl::InsecureBitGen rng_insecure_bitgen { };
    bench_bits_stl<absl::InsecureBitGen>(rng_insecure_bitgen, "absl::InsecureBitGen");
}

[[maybe_unused]] void pcg_main()
{
    // All tests were passed
    pcg32 rng_pcg32 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg32>(rng_pcg32, "PCG::pcg32");

    // All tests were passed
    pcg64 rng_pcg64 { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg64>(rng_pcg64, "PCG::pcg64");

    // All tests were passed
    pcg32_fast rng_pcg32_fast { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg32_fast>(rng_pcg32_fast, "PCG::pcg32_fast");

    // All tests were passed
    pcg64_fast rng_pcg64_fast { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg64_fast>(rng_pcg64_fast, "PCG::pcg64_fast");

    // All tests were passed
    pcg32_oneseq_once_insecure rng_pcg32_oneseq_once_insecure { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg32_oneseq_once_insecure>(rng_pcg32_oneseq_once_insecure, "PCG::pcg32_oneseq_once_insecure");

    // All tests were passed
    pcg64_oneseq_once_insecure rng_pcg64_oneseq_once_insecure { static_cast<unsigned int>(seed()) };
    bench_bits_stl<pcg64_oneseq_once_insecure>(rng_pcg64_oneseq_once_insecure, "PCG::pcg64_oneseq_once_insecure");
}

[[maybe_unused]] void mkl_main()
{
    // All tests were passed
    bench_bits_mkl(VSL_BRNG_MT19937, "MKL::VSL_BRNG_MT19937");

    // All tests were passed
    bench_bits_mkl(VSL_BRNG_MT2203, "MKL::VSL_BRNG_MT2203");

    // All tests were passed
    bench_bits_mkl(VSL_BRNG_SFMT19937, "MKL::VSL_BRNG_SFMT19937");

    // All tests were passed
    bench_bits_mkl(VSL_BRNG_ARS5, "MKL::VSL_BRNG_ARS5");

    // All tests were passed
    bench_bits_mkl(VSL_BRNG_PHILOX4X32X10, "MKL::VSL_BRNG_PHILOX4X32X10");

    // All tests were passed
    bench_bits_mkl(VSL_BRNG_NONDETERM, "MKL::VSL_BRNG_NONDETERM");

    // These generators cannot generate 32-bit integers directly.
    // bench_bits_mkl(VSL_BRNG_MCG31, "MKL::VSL_BRNG_MCG31");
    // bench_bits_mkl(VSL_BRNG_R250, "MKL::VSL_BRNG_R250");
    // bench_bits_mkl(VSL_BRNG_SOBOL, "MKL::VSL_BRNG_SOBOL");
    // bench_bits_mkl(VSL_BRNG_MRG32K3A, "MKL::VSL_BRNG_MRG32K3A");
    // bench_bits_mkl(VSL_BRNG_NIEDERR, "MKL::VSL_BRNG_NIEDERR");

    // These things are abstract base classes.
    // bench_bits_mkl(VSL_BRNG_IABSTRACT, "MKL::VSL_BRNG_IABSTRACT");
    // bench_bits_mkl(VSL_BRNG_DABSTRACT, "MKL::VSL_BRNG_DABSTRACT");
    // bench_bits_mkl(VSL_BRNG_SABSTRACT, "MKL::VSL_BRNG_SABSTRACT");

    // These two are deprecated.
    // bench_bits_mkl(VSL_BRNG_WH, "MKL::VSL_BRNG_WH");
    // bench_bits_mkl(VSL_BRNG_MCG59, "MKL::VSL_BRNG_MCG59");
}

[[maybe_unused]] void vmt19937_main()
{
    // All tests were passed
    VMT19937RandomDevice rng_vmt19937_random_device { };
    bench_bits_stl<decltype(rng_vmt19937_random_device)>(rng_vmt19937_random_device, "VMT19937RandomDevice");

    // All tests were passed
    VSFMT19937RandomDevice rng_vsfmt19937_random_device { };
    bench_bits_stl<decltype(rng_vsfmt19937_random_device)>(rng_vsfmt19937_random_device, "VSFMT19937RandomDevice");
}

[[maybe_unused]] void gsl_main()
{
    // All tests were passed
    bench_gsl(gsl_rng_mt19937);

    // All tests were passed
    bench_gsl(gsl_rng_mt19937_1999);

    // All tests were passed
    bench_gsl(gsl_rng_mt19937_1998);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_ranlxs0);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_ranlxs1);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_ranlxs2);

    // All tests were passed
    bench_gsl(gsl_rng_ranlxd1);

    // All tests were passed
    bench_gsl(gsl_rng_ranlxd2);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_ranlux);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_ranlux389);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_cmrg);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_mrg);

    // All tests were passed
    bench_gsl(gsl_rng_taus);

    // All tests were passed
    bench_gsl(gsl_rng_taus2);

    // All tests were passed
    bench_gsl(gsl_rng_gfsr4);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_rand);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_random_bsd);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_random_libc5);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_random_glibc2);

    //        Test                          p-value
    // ----------------------------------------------
    //  1  BirthdaySpacings                 eps
    //  2  Collision                      1 - eps1
    // ----------------------------------------------
    // All other tests were passed
    // bench_gsl(gsl_rng_rand48);

    //        Test                          p-value
    // ----------------------------------------------
    //  1  BirthdaySpacings                 eps
    //  2  Collision                      1 - eps1
    //  6  MaxOft AD                      3.2e-10
    // ----------------------------------------------
    // All other tests were passed
    // bench_gsl(gsl_rng_ranf);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_ranmar);

    //        Test                          p-value
    // ----------------------------------------------
    //  3  Gap                              eps
    // ----------------------------------------------
    // All other tests were passed
    // bench_gsl(gsl_rng_r250);

    //        Test                          p-value
    // ----------------------------------------------
    // 10  RandomWalk1 R                   0.9997
    // ----------------------------------------------
    // All other tests were passed
    // bench_gsl(gsl_rng_tt800);

    //        Test                          p-value
    // ----------------------------------------------
    //  1  BirthdaySpacings                 eps
    //  2  Collision                        eps
    //  6  MaxOft                           eps
    //  6  MaxOft AD                      1 - eps1
    // 10  RandomWalk1 H                    eps
    // 10  RandomWalk1 M                    eps
    // 10  RandomWalk1 J                    eps
    // 10  RandomWalk1 R                    eps
    // 10  RandomWalk1 C                    eps
    // ----------------------------------------------
    // All other tests were passed
    // bench_gsl(gsl_rng_vax);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_transputer);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_randu);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_minstd);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_uni);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_uni32);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_slatec);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_zuf);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_knuthran2);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_knuthran2002);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_knuthran);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_borosh13);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_fishman18);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_fishman20);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_lecuyer21);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_waterman14);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_fishman2x);

    // Range too small, skipping.
    // bench_gsl(gsl_rng_coveyou);
}

[[maybe_unused]] void other_rngs_main()
{
    // All tests were passed
    arc4_rand32 rng_arc4 { static_cast<arc4_rand32::result_type>(seed()) };
    bench_bits_stl<arc4_rand32>(rng_arc4, "others::arc4_rand32");

    // All tests were passed
    arc4_rand64 rng_arc4_64 { static_cast<arc4_rand64::result_type>(seed()) };
    bench_bits_stl<arc4_rand64>(rng_arc4_64, "others::arc4_rand64");

    // All tests were passed
    gjrand32 rng_gjrand32 { static_cast<gjrand32::result_type>(seed()) };
    bench_bits_stl<gjrand32>(rng_gjrand32, "others::gjrand32");

    // All tests were passed
    gjrand64 rng_gjrand63 { static_cast<gjrand64::result_type>(seed()) };
    bench_bits_stl<gjrand64>(rng_gjrand63, "others::gjrand64");

    // All tests were passed
    jsf32 rng_jsf32 { static_cast<jsf32::result_type>(seed()) };
    bench_bits_stl<jsf32>(rng_jsf32, "others::jsf32");

    //        Test                          p-value
    // ----------------------------------------------
    // 10  RandomWalk1 M                   0.9996
    // ----------------------------------------------
    // All other tests were passed
    // jsf64 rng_jsf64 { static_cast<jsf64::result_type>(seed()) };
    // bench_bits_stl<jsf64>(rng_jsf64, "others::jsf64");

    // All tests were passed
    mcg128 rng_mcg128 { static_cast<mcg128::result_type>(seed()) };
    bench_bits_stl<mcg128>(rng_mcg128, "others::mcg128");

    // All tests were passed
    mcg128_fast rng_mcg128_fast { static_cast<mcg128_fast::result_type>(seed()) };
    bench_bits_stl<mcg128_fast>(rng_mcg128_fast, "others::mcg128_fast");

    // All tests were passed
    sfc32 rng_sfc32 { static_cast<sfc32::result_type>(seed()) };
    bench_bits_stl<sfc32>(rng_sfc32, "others::sfc32");

    // All tests were passed
    sfc64 rng_sfc64 { static_cast<sfc64::result_type>(seed()) };
    bench_bits_stl<sfc64>(rng_sfc64, "others::sfc64");

    // All tests were passed
    splitmix32 rng_splitmix32 { static_cast<splitmix32::result_type>(seed()) };
    bench_bits_stl<splitmix32>(rng_splitmix32, "others::splitmix32");

    // All tests were passed
    splitmix64 rng_splitmix64 { static_cast<splitmix64::result_type>(seed()) };
    bench_bits_stl<splitmix64>(rng_splitmix64, "others::splitmix64");
}

[[maybe_unused]] void xso_main()
{
    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_2x32_star, uint32_t> x01 { };
    bench_bits_stl<decltype(x01)>(x01, "xoroshiro::2x32*");

    //        Test                          p-value
    // ----------------------------------------------
    // 10  RandomWalk1 J                   0.9998
    // ----------------------------------------------
    // All other tests were passed
    // XoroshiroWrapper<old::xoroshiro_2x32_star_star, uint32_t> x02 {};
    // bench_bits_stl<decltype(x02)>(x02, "xoroshiro::2x32**");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_4x32_plus, uint32_t, 4> x03 { };
    bench_bits_stl<decltype(x03)>(x03, "xoshiro::4x32+");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_4x32_plus_plus, uint32_t, 4> x04 { };
    bench_bits_stl<decltype(x04)>(x04, "xoshiro::4x32++");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_4x32_star_star, uint32_t, 4> x05 { };
    bench_bits_stl<decltype(x05)>(x05, "xoshiro::4x32**");

    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_2x64_plus, uint64_t, 2> x06 { };
    bench_bits_stl<decltype(x06)>(x06, "xoroshiro::2x64+");

    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_2x64_plus_plus, uint64_t, 2> x07 { };
    bench_bits_stl<decltype(x07)>(x07, "xoroshiro::2x64++");

    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_2x64_star_star, uint64_t, 2> x08 { };
    bench_bits_stl<decltype(x08)>(x08, "xoroshiro::2x64**");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_4x64_plus, uint64_t, 4> x09 { };
    bench_bits_stl<decltype(x09)>(x09, "xoshiro::4x64+");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_4x64_plus_plus, uint64_t, 4> x10 { };
    bench_bits_stl<decltype(x10)>(x10, "xoshiro::4x64++");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_4x64_star_star, uint64_t, 4> x11 { };
    bench_bits_stl<decltype(x11)>(x11, "xoshiro::4x64**");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_8x64_plus, uint64_t, 8> x12 { };
    bench_bits_stl<decltype(x12)>(x12, "xoshiro::8x64+");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_8x64_plus_plus, uint64_t, 8> x13 { };
    bench_bits_stl<decltype(x13)>(x13, "xoshiro::8x64++");

    // All tests were passed
    XoroshiroWrapper<old::xoshiro_8x64_star_star, uint64_t, 8> x14 { };
    bench_bits_stl<decltype(x14)>(x14, "xoshiro::8x64**");

    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_16x64_star, uint64_t, 16> x15 { };
    bench_bits_stl<decltype(x15)>(x15, "xoroshiro::16x64*");

    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_16x64_star_star, uint64_t, 16> x16 { };
    bench_bits_stl<decltype(x16)>(x16, "xoroshiro::16x64**");

    // All tests were passed
    XoroshiroWrapper<old::xoroshiro_16x64_plus_plus, uint64_t, 16> x17 { };
    bench_bits_stl<decltype(x17)>(x17, "xoroshiro::16x64++");
}

} // namespace

int main() noexcept
{
    // Ignore SIGPIPE to prevent crashes when the subprocess exits early
    std::signal(SIGPIPE, SIG_IGN);
    // Clear previous results
    std::ofstream { "testu01_results.txt", std::ios::binary }.close();
    // Executable testu01_main should be located in the same directory of this executable.
    testu01_path = (boost::filesystem::read_symlink("/proc/self/exe").parent_path() / "testu01_main").string();
    if (!boost::filesystem::exists(testu01_path) || !boost::filesystem::is_regular_file(testu01_path)) {
        std::cerr << "Error: testu01_main executable not found in the same directory.\n";
        return EXIT_FAILURE;
    }
    stl_main();
    absl_main();
    pcg_main();
    boost_main();
    mkl_main();
    vmt19937_main();
    gsl_main();
    xso_main();
    other_rngs_main();

    return EXIT_SUCCESS;
}
