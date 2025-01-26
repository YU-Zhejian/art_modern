#pragma once

#include "libam/bam/BamLFIO.hh"
#include "libam/bam/BamOptions.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseFileReadOutput.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <htslib/sam.h>

#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

class HeadlessBamReadOutput : public BaseFileReadOutput {
public:
    DELETE_MOVE(HeadlessBamReadOutput)
    DELETE_COPY(HeadlessBamReadOutput)

    bool require_alignment() const override;

    HeadlessBamReadOutput(const std::string& filename, const BamOptions& sam_options);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    ~HeadlessBamReadOutput() override;

private:
    samFile* sam_file_;
    sam_hdr_t* sam_header_;
    const BamOptions sam_options_;
    BamLFIO lfio_;
};

class HeadlessBamReadOutputFactory : public BaseReadOutputFactory {
public:
    DELETE_MOVE(HeadlessBamReadOutputFactory)
    DELETE_COPY(HeadlessBamReadOutputFactory)
    HeadlessBamReadOutputFactory() = default;

    [[nodiscard]] const std::string name() const override { return "HeadlessBam"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
        const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const override;
    ~HeadlessBamReadOutputFactory() override;

private:
    BamOptions sam_options_;
};

} // namespace labw::art_modern
