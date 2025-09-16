#pragma once

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/lockfree/SimpleLFIO.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>

#include <memory>
#include <string>

namespace labw::art_modern {

class FastaReadOutput : public BaseReadOutput {
public:
    DELETE_MOVE(FastaReadOutput)
    DELETE_COPY(FastaReadOutput)
    explicit FastaReadOutput(const std::string& filename, int n_threads);
    void writeSE(const ProducerToken& token, const PairwiseAlignment& pwa) override;
    void writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    bool require_alignment() const override;
    ProducerToken get_producer_token() override;
    ~FastaReadOutput() override;

private:
    SimpleLFIO lfio_;
};

class FastaReadOutputFactory : public BaseReadOutputFactory {
public:
    DELETE_MOVE(FastaReadOutputFactory)
    DELETE_COPY(FastaReadOutputFactory)
    FastaReadOutputFactory() = default;
    ~FastaReadOutputFactory() override = default;

    [[nodiscard]] const std::string name() const override { return "FASTA"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    [[nodiscard]] std::shared_ptr<BaseReadOutput> create(const OutParams& params) const override;
};
} // namespace labw::art_modern
