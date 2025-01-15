#include "fasta/BaseFastaFetch.hh"
#include "fasta/FaidxFetch.hh"
#include "fasta/InMemoryFastaFetch.hh"

#include <htslib/hts.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>


using namespace labw::art_modern;

void bench_ff(BaseFastaFetch* ff, const std::string& name)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::tuple<int, hts_pos_t, hts_pos_t>> queries;
    hts_pos_t gen_start = 0;
    hts_pos_t gen_end = 0;

    for (int i = 0; i < ff->num_seqs(); i++) {
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

int main()
{
    std::vector<std::string> data { "ce11.mRNA_head", "ce11_chr1" };

    for (const auto& datum : data) {
        BaseFastaFetch* ff = new FaidxFetch(datum + ".fa");
        bench_ff(ff, datum + "-faidx-fa");
        delete ff;

        ff = new FaidxFetch(datum + ".fa.gz");
        bench_ff(ff, datum + "-faidx-fa.gz");
        delete ff;

        ff = new InMemoryFastaFetch(datum + ".fa");
        bench_ff(ff, datum + "-memory-fa");
        delete ff;
    }
}
