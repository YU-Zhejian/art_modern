#pragma once

#include "libam/ds/PairwiseAlignment.hh"
#include "libam/lockfree/SimpleLFIO.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

class PwaReadOutput : public BaseReadOutput {
public:
    DELETE_MOVE(PwaReadOutput)
    DELETE_COPY(PwaReadOutput)

    explicit PwaReadOutput(const std::string& filename, const std::vector<std::string>& args);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

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
    std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
        const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const override;
};
} // namespace labw::art_modern
