#include "libam_support/Constants.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/FastqReadOutput.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/si_utils.hh"

#include "align_blkring.hh"

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace labw::art_modern;

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
constexpr int NTHREAD = 20;

void working_thread(const std::shared_ptr<BaseReadOutput>& bro, const std::size_t num_records)
{
    const auto token = bro->get_producer_token();
    for (std::size_t i = 0; i < num_records; i++) {
        bro->writeSE(token, pwa);
    }
}

void bench(const std::shared_ptr<BaseReadOutput>& bro, const std::string& name, const std::size_t nthread)
{
    std::cout << "Benchmarking " << name << " with " << nthread << " threads" << std::endl;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < nthread; i++) {
        std::thread t(working_thread, bro, (200ULL * M_SIZE) / nthread);
        threads.emplace_back(std::move(t));
    }
    for (auto& t : threads) {
        t.join();
    }
    bro->close();
    end = std::chrono::high_resolution_clock::now();
}

void speed2devnull()
{
    const std::size_t block_size = 4ULL * K_SIZE;
    const std::size_t n_blocks = 20ULL * M_SIZE;
    const std::string block = std::string(block_size, '\0');
    auto start = std::chrono::high_resolution_clock::now();
    auto ofs = std::ofstream(DEVNULL, std::ios::out | std::ios::binary);
    for (std::size_t i = 0; i < n_blocks; i++) {
        ofs << block;
    }
    ofs.close();
    auto end = std::chrono::high_resolution_clock::now();
    const double speed = 1.0 * (n_blocks * block_size)
        / std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * 1000;
    std::cout << "Speed: " << to_si(speed) << "B/s" << std::endl;
}

} // namespace

int main()
{
    speed2devnull();
    auto iss = std::istringstream { fasta };
    const auto ff = std::make_shared<InMemoryFastaFetch>(iss);

    // Temporarily disable logging
    bench(std::make_shared<FastqReadOutput>(DEVNULL, NTHREAD), "FastqReadOutput", NTHREAD);

    return EXIT_SUCCESS;
}
