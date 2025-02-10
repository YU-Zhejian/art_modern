#include "libam_support/Constants.hh"
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/out/BamReadOutput.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/FastaReadOutput.hh"
#include "libam_support/out/FastqReadOutput.hh"
#include "libam_support/out/HeadlessBamReadOutput.hh"
#include "libam_support/out/PwaReadOutput.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <atomic>
#include <fmt/core.h>

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace logging = boost::log;
using namespace labw::art_modern;

#if 0
template <typename T> class LockedIO {
public:
    DELETE_MOVE(LockedIO)
    DELETE_COPY(LockedIO)

    constexpr static const std::chrono::duration sleep_time = std::chrono::microseconds(10);

    explicit LockedIO(std::string name)
        : name_(std::move(name)) { };

    virtual ~LockedIO() = default;

    void push(T&& value)
    {
        std::scoped_lock lock(mutex_);
        num_reads_in_++;
        write(std::move(value));
        num_reads_out_++;
    }
    void start() { start_time_ = std::chrono::high_resolution_clock::now(); }
    virtual void flush_and_close() { };

    void stop()
    {
        flush_and_close();
        end_time_ = std::chrono::high_resolution_clock::now();
        if (!had_logged_) {
            log_();
        }
        had_logged_ = true;
    }

    virtual void write(T value) = 0;

protected:
    std::atomic<std::size_t> num_bytes_out_ = 0;
    const std::string name_;

private:
    std::chrono::high_resolution_clock::time_point start_time_;
    std::chrono::high_resolution_clock::time_point end_time_;
    std::atomic<std::size_t> num_reads_in_ = 0;
    std::atomic<std::size_t> num_reads_out_ = 0;
    std::atomic<std::size_t> num_wait_in_ = 0;
    std::atomic<std::size_t> num_wait_out_not_full_ = 0;
    std::atomic<std::size_t> num_wait_out_empty_ = 0;
    std::atomic<std::size_t> num_nowait_out_ = 0;
    std::atomic<bool> had_logged_ = false;

    std::mutex mutex_;

    void run()
    {
        // Does nothing!
    }
    void log_() const
    {
        const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_).count();
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: Finished, consuming " << num_reads_in_ << " reads and writes "
                                << num_reads_out_ << " reads.";
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: N. Waitings (I/ONotFull/OEmpty): " << num_wait_in_ << " / "
                                << num_wait_out_not_full_ << "("
                                << (100.0 * num_wait_out_not_full_
                                       / (num_wait_out_empty_ + num_wait_out_not_full_ + num_nowait_out_))
                                << "%) / " << num_wait_out_empty_ << "("
                                << (100.0 * num_wait_out_empty_
                                       / (num_wait_out_empty_ + num_wait_out_not_full_ + num_nowait_out_))
                                << "%).";
        BOOST_LOG_TRIVIAL(info) << name_ << " LockFreeIO: " << to_si(num_bytes_out_) << "B written in " << time / 1000.0
                                << " seconds. Speed: " << to_si(1.0 * num_bytes_out_ / (time / 1000.0)) << "B/s.";
    }
};
#endif

class EmptyLFIO : public LockFreeIO<std::unique_ptr<std::nullptr_t>> {
public:
    DELETE_COPY(EmptyLFIO)
    DELETE_MOVE(EmptyLFIO)
    EmptyLFIO()
        : LockFreeIO<std::unique_ptr<std::nullptr_t>>("Empty")
    {
    }
    ~EmptyLFIO() override { stop(); };
    void write([[maybe_unused]] std::unique_ptr<std::nullptr_t> /** value **/) override
    {
        // Do nothing!
    }
};

/**
 * @brief An output that passes `std::unique_ptr<std::nullptr_t>` to an LFIO.
 *
 * Designed to test the overhead caused by Moody Camel Queue and the LockFreeIO.
 */
class EmptyLFIOReadOutput final : public BaseReadOutput {
public:
    DELETE_COPY(EmptyLFIOReadOutput)
    DELETE_MOVE(EmptyLFIOReadOutput)

    EmptyLFIOReadOutput(const int nthreads) {
        lfio_.init_queue(nthreads, 0); lfio_.start();
    }
    void writeSE(const moodycamel::ProducerToken& token, [[maybe_unused]] const PairwiseAlignment& /** pwa **/) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>(), token);
    }
    void writePE(const moodycamel::ProducerToken& token, [[maybe_unused]] const PairwiseAlignment&,
        [[maybe_unused]] const PairwiseAlignment&) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>(), token);
        lfio_.push(std::make_unique<std::nullptr_t>(), token);
    }
    void close() override
    {
        lfio_.flush_and_close();
        lfio_.stop();
    }

    bool require_alignment() const override { return false; }

    ~EmptyLFIOReadOutput() override { close(); }
    moodycamel::ProducerToken get_producer_token() override { return lfio_.get_producer_token(); }

private:
    EmptyLFIO lfio_;
};
class EmptyImplicitLFIOReadOutput final : public BaseReadOutput {
public:
    DELETE_COPY(EmptyImplicitLFIOReadOutput)
    DELETE_MOVE(EmptyImplicitLFIOReadOutput)

    EmptyImplicitLFIOReadOutput(const int nthreads) {
        lfio_.init_queue(0, nthreads); lfio_.start();
    }
    void writeSE([[maybe_unused]] const moodycamel::ProducerToken& token, [[maybe_unused]] const PairwiseAlignment& /** pwa **/) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>());
    }
    void writePE([[maybe_unused]] const moodycamel::ProducerToken& token, [[maybe_unused]] const PairwiseAlignment&,
                 [[maybe_unused]] const PairwiseAlignment&) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>());
        lfio_.push(std::make_unique<std::nullptr_t>());
    }
    void close() override
    {
        lfio_.flush_and_close();
        lfio_.stop();
    }

    bool require_alignment() const override { return false; }

    ~EmptyImplicitLFIOReadOutput() override { close(); }
    moodycamel::ProducerToken get_producer_token() override { return lfio_.get_producer_token(); }

private:
    EmptyLFIO lfio_;
};


namespace {
const std::string DEVNULL = "/dev/null";
const std::string fasta = ">chr1\nGGGCGTGTTCCTGTCGGGTAACACCACCATAGCAAAGCGATTGTTTATTTGACGAGTAAGGGAGGTCATTTCTATGACGGGGGGA"
                          "CCAGAGCCGCGGTGCATCACTCTAGAACTCCAGCTTATTTACAACATGGTGAGATGATTAGATGG";
const PairwiseAlignment pwa { "read_1", "chr1",
    "GGGCGTGTTCCTGTCGGGTAACACCACCATAGCAAAGCGATTGTTTATTTGACGAGTAAG"
    "GGAGGTCATTTCTATGACGGGGGGACCAGAGCCGCGGTGCATCACTCTAGAACT"
    "CCAGCTTATTTACAACATGGTGAGATGATTAGATGG",
    "GGGCGTGTTCCTGTCGGGTAACACCACCATAGCAAGCGATTGTTTATTTGACGAGTAAGGG"
    "AAAAGGTCATTTCCATGACGGGGGGACCAGAGCCGCGGTGCATCACTCTAGAG"
    "CTCCAGCTTATTTACAACATGGTGAGATTTAGATGG",
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
    "GGGCGTGTTCCTGTCGGGTAACACCACCATAGCAAAGCGATTGTTTATTTGACGAGTAAG"
    "GG---AGGTCATTTCTATGACGGGGGGACCAGAGCCGCGGTGCATCACTCTAGAACTCCA"
    "GCTTATTTACAACATGGTGAGATGATTAGATGG",
    "GGGCGTGTTCCTGTCGGGTAACACCACCATAGC-AAGCGATTGTTTATTTGACGAGTAAG"
    "GGAAAAGGTCATTTCCATGACGGGGGGACCAGAGCCGCGGTGCATCACTCTAGAGCTCCA"
    "GCTTATTTACAACATGGTGAGAT--TTAGATGG",
    0, true };

void working_thread(const std::shared_ptr<BaseReadOutput>& bro, const std::size_t num_records)
{
    const auto token = bro->get_producer_token();
    for (std::size_t i = 0; i < num_records; i++) {
        bro->writeSE(token, pwa);
    }
}

void bench(const std::shared_ptr<BaseReadOutput>& bro, const std::string& name, std::size_t nthread, std::ostream& oss)
{
    std::cout << "Benchmarking " << name << " with " << nthread << " threads" << std::endl;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < nthread; i++) {
        std::thread t(working_thread, bro, (M_SIZE) / nthread);
        threads.emplace_back(std::move(t));
    }
    for (auto& t : threads) {
        t.join();
    }
    bro->close();
    end = std::chrono::high_resolution_clock::now();
    oss << fmt::format(
        "{}\t{}\t{}\n", name, nthread, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::flush(oss);
}

std::string get_bo_name(const std::string& name, const BamOptions& bo)
{
    return fmt::format("{}_l={}_t={}", name, bo.compress_level, bo.hts_io_threads);
}

} // namespace

int main()
{
    auto oss = std::ofstream { "time_complexity.tsv" };
    oss << "name\tthreads\ttime_complexity" << std::endl;
    auto iss = std::istringstream { fasta };
    const auto ff = std::make_shared<InMemoryFastaFetch>(iss);
    BamOptions so;
    BamOptions const bo_defaults;
    so.write_bam = false;
    std::vector<BamOptions> bo_t;
    std::vector<BamOptions> bo_l;

    for (auto& t : std::vector<int> { 1, 2, 4, 8, 16, 32, 64 }) {
        bo_t.emplace_back();
        bo_t.back().hts_io_threads = t;
    }
    for (auto& t : std::vector<char> { '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 'u' }) {
        bo_l.emplace_back();
        bo_l.back().compress_level = t;
    }
    // Temporarily disable logging
    logging::core::get()->set_filter(logging::trivial::severity > logging::trivial::fatal);
    for (int i = 0; i < 5; i++) {
        for (std::size_t const nthread : std::vector<std::size_t> { 1, 2, 4, 8, 16, 32, 64 }) {
            for (auto& bo : bo_l) {
                bench(std::make_shared<BamReadOutput>(DEVNULL, ff, bo, nthread), get_bo_name("BamReadOutput", bo),
                    nthread, oss);
            }
            for (auto& bo : bo_t) {
                bench(std::make_shared<BamReadOutput>(DEVNULL, ff, bo, nthread), get_bo_name("BamReadOutput", bo),
                    nthread, oss);
            }
            bench(std::make_shared<HeadlessBamReadOutput>(DEVNULL, bo_defaults, nthread),
                get_bo_name("HeadlessBamReadOutput", bo_defaults), nthread, oss);
            bench(std::make_shared<BamReadOutput>(DEVNULL, ff, so, nthread), get_bo_name("SamReadOutput", so), nthread,
                oss);
            bench(std::make_shared<FastqReadOutput>(DEVNULL, nthread), "FastqReadOutput", nthread, oss);
            bench(std::make_shared<FastaReadOutput>(DEVNULL, nthread), "FastaReadOutput", nthread, oss);
            bench(std::make_shared<EmptyLFIOReadOutput>(nthread), "EmptyLFIOReadOutput", nthread, oss);
            bench(std::make_shared<EmptyImplicitLFIOReadOutput>(nthread), "EmptyImplicitLFIOReadOutput", nthread, oss);
            bench(std::make_shared<PwaReadOutput>(DEVNULL, std::vector<std::string> {}, nthread), "PwaReadOutput",
                nthread, oss);
        }
    }
    logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);

    return EXIT_SUCCESS;
}
