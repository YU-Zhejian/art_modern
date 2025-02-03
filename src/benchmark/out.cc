#include "libam/Constants.hh"
#include "libam/bam/BamOptions.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/lockfree/LockFreeIO.hh"
#include "libam/out/BamReadOutput.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/DumbReadOutput.hh"
#include "libam/out/FastaReadOutput.hh"
#include "libam/out/FastqReadOutput.hh"
#include "libam/out/HeadlessBamReadOutput.hh"
#include "libam/out/PwaReadOutput.hh"
#include "libam/ref/fetch/InMemoryFastaFetch.hh"

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#include <chrono>
#include <cstdlib>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace logging = boost::log;
using namespace labw::art_modern;


class EmptyLFIOReadOutput final : public BaseReadOutput {
public:
    DELETE_COPY(EmptyLFIOReadOutput)
    DELETE_MOVE(EmptyLFIOReadOutput)

    EmptyLFIOReadOutput() {
        lfio_.start();
    }
    void writeSE( [[maybe_unused]] const PairwiseAlignment& /** pwa **/) override{
        lfio_.push(std::make_unique<nullptr_t>());
    }
    void writePE([[maybe_unused]] const PairwiseAlignment& ,[[maybe_unused]]  const PairwiseAlignment& ) override{
        lfio_.push(std::make_unique<nullptr_t>());
        lfio_.push(std::make_unique<nullptr_t>());
    }
    void close() override{
        lfio_.flush_and_close();
        lfio_.stop();
    }

    bool require_alignment() const override{
        return false;
    }

    ~EmptyLFIOReadOutput() override{
        close();
    }
private:
    class EmptyLFIO: public LockFreeIO<std::unique_ptr<nullptr_t> >{
    public:
        DELETE_COPY(EmptyLFIO)
        DELETE_MOVE(EmptyLFIO)
        EmptyLFIO() : LockFreeIO<std::unique_ptr<nullptr_t>>("Empty") {}
        void write(std::unique_ptr<nullptr_t> value) override {
            // Do nothing!
        }
    };
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
std::unordered_map<std::string, int> time_complexity {};

void bench(std::unique_ptr<BaseReadOutput> bro, const std::string& name)
{
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < M_SIZE; i++) {
        bro->writeSE(pwa);
    }
    bro->close();
    end = std::chrono::high_resolution_clock::now();
    time_complexity.emplace(name, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
}

} // namespace

int main()
{
    auto iss = std::istringstream { fasta };
    const auto ff = InMemoryFastaFetch(iss);
    BamOptions bo9;
    BamOptions bo1;
    BamOptions bou;
    BamOptions so;
    so.write_bam = false;
    bo9.compress_level = '9';
    bo1.compress_level = '1';
    bou.compress_level = 'u';
    // Temporarily disable logging
    logging::core::get()->set_filter(logging::trivial::severity > logging::trivial::fatal);

    bench(std::make_unique<BamReadOutput>(DEVNULL, &ff, bo1), "BamReadOutput (l=1)");
    bench(std::make_unique<BamReadOutput>(DEVNULL, &ff, bo9), "BamReadOutput (l=9)");
    bench(std::make_unique<BamReadOutput>(DEVNULL, &ff, bou), "BamReadOutput (l=u)");
    bench(std::make_unique<HeadlessBamReadOutput>(DEVNULL, bo1), "HeadlessBamReadOutput (l=1)");
    bench(std::make_unique<BamReadOutput>(DEVNULL, &ff, so), "SamReadOutput");
    bench(std::make_unique<FastqReadOutput>(DEVNULL), "FastqReadOutput");
    bench(std::make_unique<FastaReadOutput>(DEVNULL), "FastaReadOutput");
    bench(std::make_unique<DumbReadOutput>(), "DumbReadOutput");
    bench(std::make_unique<EmptyLFIOReadOutput>(), "EmptyLFIOReadOutput");
    bench(std::make_unique<PwaReadOutput>(DEVNULL, std::vector<std::string> {}), "PwaReadOutput");

    logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);
    for (const auto& [name, time] : time_complexity) {
        BOOST_LOG_TRIVIAL(info) << name << ": " << time << " ms";
    }
    return EXIT_SUCCESS;
}
