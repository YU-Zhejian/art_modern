#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include "ds/PairwiseAlignment.hh"
#include "lockfree/SimpleLFIO.hh"
#include "out/BaseFileReadOutput.hh"
#include "out/BaseReadOutput.hh"
#include "ref/fetch/BaseFastaFetch.hh"

namespace labw::art_modern {

class PwaReadOutput : public BaseFileReadOutput {
public:
    PwaReadOutput(PwaReadOutput&& other) = delete;
    PwaReadOutput(const PwaReadOutput&) = delete;
    PwaReadOutput& operator=(PwaReadOutput&&) = delete;
    PwaReadOutput& operator=(const PwaReadOutput&) = delete;

    explicit PwaReadOutput(const std::string& filename, const std::vector<std::string>& args);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

    void close() override;

    ~PwaReadOutput() override;

private:
    std::ofstream file_;
    SimpleLFIO lfio_;
};

class PwaReadOutputFactory : public BaseReadOutputFactory {
public:
    const std::string name() const override { return "PWA"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
};
} // namespace labw::art_modern
