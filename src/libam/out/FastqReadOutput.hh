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

class FastqReadOutput : public BaseReadOutput {
public:
    DELETE_MOVE(FastqReadOutput)
    DELETE_COPY(FastqReadOutput)
    explicit FastqReadOutput(const std::string& filename);
    void writeSE(const PairwiseAlignment& pwa) override;
    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    void close() override;
    ~FastqReadOutput() override;

    bool require_alignment() const override;

private:
    SimpleLFIO lfio_;
};

class FastqReadOutputFactory : public BaseReadOutputFactory {
public:
    DELETE_MOVE(FastqReadOutputFactory)
    DELETE_COPY(FastqReadOutputFactory)
    FastqReadOutputFactory() = default;
    ~FastqReadOutputFactory() override = default;

    [[nodiscard]] const std::string name() const override { return "FASTQ"; }
    void patch_options(boost::program_options::options_description& desc) const override;
    std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
        const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const override;
};
} // namespace labw::art_modern
