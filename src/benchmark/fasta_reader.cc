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

#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/ref/fetch/FaidxFetch.hh"
#include "libam_support/ref/fetch/InMemoryFastaFetch.hh"

#include <htslib/hts.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace labw::art_modern; // NOLINT
namespace {
void bench_ff(const std::unique_ptr<BaseFastaFetch>& ff, const std::string& name)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::tuple<int, hts_pos_t, hts_pos_t>> queries;
    hts_pos_t gen_start = 0;
    hts_pos_t gen_end = 0;

    for (std::size_t i = 0; i < ff->num_seqs(); i++) {
        std::uniform_int_distribution<hts_pos_t> coord_dist(0, ff->seq_len(i));
        for (int j = 0; j < 1000; j++) {
            gen_start = coord_dist(gen);
            gen_end = coord_dist(gen);
            if (gen_start == gen_end) {
                continue;
            }
            queries.emplace_back(i, std::min(gen_start, gen_end), std::max(gen_start, gen_end));
        }
    }

    size_t fetched_size = 0;
    for (const auto& [contig_id, start, end] : queries) {
        fetched_size += ff->fetch(contig_id, start, end).size();
    }
    std::cout << name << " fetched_size: " << fetched_size << std::endl;
}
} // namespace
int main()
{
    const std::vector<std::string> data { "ce11.mRNA_head", "ce11_chr1" };

    for (const auto& datum : data) {
        std::unique_ptr<BaseFastaFetch> ff = std::make_unique<FaidxFetch>(datum + ".fa");
        bench_ff(std::move(ff), datum + "-faidx-fa");

        ff = std::make_unique<FaidxFetch>(datum + ".fa.gz");
        bench_ff(std::move(ff), datum + "-faidx-fa.gz");

        ff = std::make_unique<InMemoryFastaFetch>(datum + ".fa");
        bench_ff(std::move(ff), datum + "-memory-fa");
    }
}
