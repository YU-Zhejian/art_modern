#pragma once

#include <fstream>
#include <iostream>
#include <mutex>

#include "out/BaseReadOutput.hh"

namespace labw::art_modern {
class PwaReadOutput : public BaseReadOutput {
public:
    PwaReadOutput(PwaReadOutput&& other) = delete;
    PwaReadOutput(const PwaReadOutput&) = delete;
    PwaReadOutput& operator=(PwaReadOutput&&) = delete;
    PwaReadOutput& operator=(const PwaReadOutput&) = delete;

    explicit PwaReadOutput(const std::string& filename);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

    void close() override;

    ~PwaReadOutput() override;

private:
    std::ofstream file_;
    std::mutex mutex_;
    const std::string filename;
    bool is_closed_;
};

class PwaReadOutputFactory : public BaseReadOutputFactory {
public:
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
};
}
