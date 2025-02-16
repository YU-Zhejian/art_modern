#pragma once

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/SimpleLFIO.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <concurrentqueue.h>

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
    explicit FastqReadOutput(const std::string& filename, int n_threads);
    void writeSE(const moodycamel::ProducerToken& token, const PairwiseAlignment& pwa) override;
    void writePE(
        const moodycamel::ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;
    moodycamel::ProducerToken get_producer_token() override;
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
    std::shared_ptr<BaseReadOutput> create(const OutParams& params) const override;
};
} // namespace labw::art_modern
