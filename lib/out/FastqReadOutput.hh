#pragma once

#include <fstream>
#include <mutex>

#include "out/BaseReadOutput.hh"

namespace labw::art_modern {
class FastqReadOutput : public BaseReadOutput {
public:
    FastqReadOutput(FastqReadOutput&& other) = delete;
    FastqReadOutput(const FastqReadOutput&) = delete;
    FastqReadOutput& operator=(FastqReadOutput&&) = delete;
    FastqReadOutput& operator=(const FastqReadOutput&) = delete;

    explicit FastqReadOutput(const std::string& filename);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

    void close() override;

    ~FastqReadOutput() override;

private:
    std::ofstream file_;
    std::mutex mutex_;
    const std::string filename;
    bool is_closed_;
};

class FastqReadOutputFactory : public BaseReadOutputFactory {
public:
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
};
}
