/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "libam_support/Constants.hh"
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/LockFreeIO.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BamReadOutput.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/FastaReadOutput.hh"
#include "libam_support/out/FastqReadOutput.hh"
#include "libam_support/out/HeadlessBamReadOutput.hh"
#include "libam_support/out/PwaReadOutput.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <fmt/format.h>

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace labw::art_modern;

class EmptyLFIO final : public LockFreeIO<std::unique_ptr<std::nullptr_t>> {
public:
    DELETE_COPY(EmptyLFIO)
    DELETE_MOVE(EmptyLFIO)
    EmptyLFIO()
        : LockFreeIO("Empty")
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

    explicit EmptyLFIOReadOutput(const int nthreads)
    {
        lfio_.init_queue(nthreads, 0);
        lfio_.start();
    }
    void writeSE(const ProducerToken& token, [[maybe_unused]] const PairwiseAlignment& /** pwa **/) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>(), token);
    }
    void writePE(const ProducerToken& token, [[maybe_unused]] const PairwiseAlignment& /** pwa1 **/,
        [[maybe_unused]] const PairwiseAlignment& /** pwa2 **/) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>(), token);
        lfio_.push(std::make_unique<std::nullptr_t>(), token);
    }
    void close() override { lfio_.stop(); }

    [[nodiscard]] bool require_alignment() const override { return false; }

    ~EmptyLFIOReadOutput() override { close(); }
    ProducerToken get_producer_token() override { return lfio_.get_producer_token(); }

private:
    EmptyLFIO lfio_;
};
class EmptyImplicitLFIOReadOutput final : public BaseReadOutput {
public:
    DELETE_COPY(EmptyImplicitLFIOReadOutput)
    DELETE_MOVE(EmptyImplicitLFIOReadOutput)

    explicit EmptyImplicitLFIOReadOutput(const int nthreads)
    {
        lfio_.init_queue(0, nthreads);
        lfio_.start();
    }
    void writeSE([[maybe_unused]] const ProducerToken& /** token **/,
        [[maybe_unused]] const PairwiseAlignment& /** pwa **/) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>());
    }
    void writePE([[maybe_unused]] const ProducerToken& /** token **/,
        [[maybe_unused]] const PairwiseAlignment& /** pwa1 **/,
        [[maybe_unused]] const PairwiseAlignment& /** pwa2 **/) override
    {
        lfio_.push(std::make_unique<std::nullptr_t>());
        lfio_.push(std::make_unique<std::nullptr_t>());
    }
    void close() override { lfio_.stop(); }

    [[nodiscard]] bool require_alignment() const override { return false; }

    ~EmptyImplicitLFIOReadOutput() override { close(); }
    ProducerToken get_producer_token() override { return lfio_.get_producer_token(); }

private:
    EmptyLFIO lfio_;
};

namespace {
const std::string DEVNULL = "/dev/null";
const std::string fasta = ">chr1\nGGGCGTGTTCCTGTCGGGTAACACCACCATAGCAAAGCGATTGTTTATTTGACGAGTAAGGGAGGTCATTTCTATGACGGGGGGA"
                          "CCAGAGCCGCGGTGCATCACTCTAGAACTCCAGCTTATTTACAACATGGTGAGATGATTAGATGG";
const std::vector<am_qual_t> QUALS(150, 0);

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
    QUALS,
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
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < nthread; i++) {
        std::thread t(working_thread, bro, M_SIZE / nthread);
        threads.emplace_back(std::move(t));
    }
    for (auto& t : threads) {
        t.join();
    }
    bro->close();
    auto end = std::chrono::high_resolution_clock::now();
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

    for (auto& t : std::vector { 1, 2, 4, 8, 16, 32, 64 }) {
        bo_t.emplace_back();
        bo_t.back().hts_io_threads = t;
    }
    for (auto& t : std::vector { '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 'u' }) {
        bo_l.emplace_back();
        bo_l.back().compress_level = t;
    }
    // Temporarily disable logging
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

    return EXIT_SUCCESS;
}
