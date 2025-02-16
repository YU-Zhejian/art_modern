#pragma once

#include "libam_support/bam/BamLFIO.hh"
#include "libam_support/bam/BamOptions.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <concurrentqueue.h>

#include <boost/program_options/options_description.hpp>

#include <htslib/sam.h>

#include <memory>
#include <string>

namespace labw::art_modern {

class BamReadOutput : public BaseReadOutput {
public:
    DELETE_MOVE(BamReadOutput)
    DELETE_COPY(BamReadOutput)

    BamReadOutput(const std::string& filename, const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
        const BamOptions& sam_options, int n_threads);
    void writeSE(const moodycamel::ProducerToken& token, const PairwiseAlignment& pwa) override;
    void writePE(
        const moodycamel::ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    ~BamReadOutput() override;

    [[nodiscard]] bool require_alignment() const override;
    moodycamel::ProducerToken get_producer_token() override;

private:
    samFile* sam_file_;
    sam_hdr_t* sam_header_;
    const BamOptions sam_options_;
    BamLFIO lfio_;
};

class BamReadOutputFactory : public BaseReadOutputFactory {
public:
    DELETE_MOVE(BamReadOutputFactory)
    DELETE_COPY(BamReadOutputFactory)
    BamReadOutputFactory() = default;

    [[nodiscard]] const std::string name() const override;
    void patch_options(boost::program_options::options_description& desc) const override;
    [[nodiscard]] std::shared_ptr<BaseReadOutput> create(const OutParams& params) const override;
    ~BamReadOutputFactory() override;

private:
    BamOptions sam_options_;
};
} // namespace labw::art_modern
