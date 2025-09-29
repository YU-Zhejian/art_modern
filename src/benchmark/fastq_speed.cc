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
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/FastqReadOutput.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"
#include "libam_support/utils/si_utils.hh"

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
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> threads;
    for (std::size_t i = 0; i < nthread; i++) {
        std::thread t(working_thread, bro, (200ULL * M_SIZE) / nthread);
        threads.emplace_back(std::move(t));
    }
    for (auto& t : threads) {
        t.join();
    }
    bro->close();
    auto end = std::chrono::high_resolution_clock::now();
}

void speed2devnull()
{
    constexpr std::size_t block_size = 4ULL * K_SIZE;
    constexpr std::size_t n_blocks = 20ULL * M_SIZE;
    const auto block = std::string(block_size, '\0');
    const auto start = std::chrono::high_resolution_clock::now();
    auto ofs = std::ofstream(DEVNULL, std::ios::out | std::ios::binary);
    for (std::size_t i = 0; i < n_blocks; i++) {
        ofs << block;
    }
    ofs.close();
    const auto end = std::chrono::high_resolution_clock::now();
    const double speed = static_cast<double>(n_blocks * block_size)
        / static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) * 1000.0;
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
