#pragma once

#include <fstream>
#include <mutex>

#include "LockFreeIO.hh"
#include "out/BaseReadOutput.hh"

namespace labw::art_modern {

class FastqLFIO : public LockFreeIO<std::stringstream> {
public:
    void write(std::stringstream* ss) override
    {
        out_ << ss->str();
        delete ss;
    }
    FastqLFIO(std::ostream& out)
        : out_(out)
    {
    }

private:
    std::ostream& out_;
};

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
    const std::string filename;
    bool is_closed_;
    FastqLFIO lfio_;
};

class FastqReadOutputFactory : public BaseReadOutputFactory {
public:
    const std::string name() const override { return "FASTQ"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
};
}
