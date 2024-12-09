#pragma once

#include <fstream>
#include <iostream>

#include "LockFreeIO.hh"
#include "out/BaseFileReadOutput.hh"

namespace labw::art_modern {
class PwaLFIO : public LockFreeIO<std::ostringstream> {
public:
    void write(std::unique_ptr<std::ostringstream> ss) override { out_ << ss->str(); }
    explicit PwaLFIO(std::ostream& out)
        : out_(out)
    {
    }

private:
    std::ostream& out_;
};

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
    PwaLFIO lfio_;
};

class PwaReadOutputFactory : public BaseReadOutputFactory {
public:
    const std::string name() const override { return "PWA"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
};
}
