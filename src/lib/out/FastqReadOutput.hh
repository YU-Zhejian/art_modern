#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "BaseReadOutput.hh"
#include "ds/PairwiseAlignment.hh"
#include "lockfree/SimpleLFIO.hh"
#include "out/BaseFileReadOutput.hh"
#include "out/BaseReadOutput.hh"
#include "ref/fetch/BaseFastaFetch.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

namespace labw::art_modern {

class FastqReadOutput : public BaseFileReadOutput {
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
    SimpleLFIO lfio_;
};

class FastqReadOutputFactory : public BaseReadOutputFactory {
public:
    [[nodiscard]] const std::string name() const override { return "FASTQ"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
};
} // namespace labw::art_modern
