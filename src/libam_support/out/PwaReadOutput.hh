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
#include <vector>

namespace labw::art_modern {

class PwaReadOutput : public BaseReadOutput {
public:
    DELETE_MOVE(PwaReadOutput)
    DELETE_COPY(PwaReadOutput)

    explicit PwaReadOutput(const std::string& filename, const std::vector<std::string>& args, int n_threads);
    void writeSE(const ProducerToken& token, const PairwiseAlignment& pwa) override;
    void writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    ProducerToken get_producer_token() override;

    bool require_alignment() const override;

    void close() override;

    ~PwaReadOutput() override;

private:
    SimpleLFIO lfio_;
};

class PwaReadOutputFactory : public BaseReadOutputFactory {
public:
    DELETE_MOVE(PwaReadOutputFactory)
    DELETE_COPY(PwaReadOutputFactory)
    PwaReadOutputFactory() = default;
    ~PwaReadOutputFactory() override = default;

    [[nodiscard]] const std::string name() const override { return "PWA"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    std::shared_ptr<BaseReadOutput> create(const OutParams& params) const override;
};
} // namespace labw::art_modern
