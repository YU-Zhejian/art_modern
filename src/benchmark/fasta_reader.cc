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
void bench_ff(std::unique_ptr<BaseFastaFetch> ff, const std::string& name)
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
