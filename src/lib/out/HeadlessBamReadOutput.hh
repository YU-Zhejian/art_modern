#pragma once

#include "BamUtils.hh"
#include "out/BaseFileReadOutput.hh"

namespace labw::art_modern {

class HeadlessBamReadOutput : public BaseFileReadOutput {
public:
    HeadlessBamReadOutput(HeadlessBamReadOutput&& other) = delete;
    HeadlessBamReadOutput(const HeadlessBamReadOutput&) = delete;
    HeadlessBamReadOutput& operator=(HeadlessBamReadOutput&&) = delete;
    HeadlessBamReadOutput& operator=(const HeadlessBamReadOutput&) = delete;

    HeadlessBamReadOutput(const std::string& filename, const SamOptions& sam_options);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    ~HeadlessBamReadOutput() override;

private:
    samFile* sam_file_;
    sam_hdr_t* sam_header_;
    const SamOptions sam_options_;
    BamLFIO lfio_;
};

class HeadlessBamReadOutputFactory : public BaseReadOutputFactory {
public:
    const std::string name() const override { return "HeadlessBam"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
    ~HeadlessBamReadOutputFactory() override;

private:
    SamOptions sam_options_;
};

}
