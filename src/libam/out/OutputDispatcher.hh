#pragma once

#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

class OutputDispatcher : public BaseReadOutput {

public:
    OutputDispatcher(OutputDispatcher&& other) = delete;
    OutputDispatcher(const OutputDispatcher&) = delete;
    OutputDispatcher& operator=(OutputDispatcher&&) = delete;
    OutputDispatcher& operator=(const OutputDispatcher&) = delete;

    OutputDispatcher() = default;

    void add(BaseReadOutput* output);

    void writeSE(const PairwiseAlignment& pwa) override;

    void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) override;

    void close() override;

    ~OutputDispatcher() override;

private:
    std::vector<BaseReadOutput*> outputs_;
};

class OutputDispatcherFactory : public BaseReadOutputFactory {
public:
    DELETE_MOVE(OutputDispatcherFactory)
    DELETE_COPY(OutputDispatcherFactory)

    OutputDispatcherFactory();

    [[nodiscard]] const std::string name() const override { return "OD"; }
    void add(std::shared_ptr<BaseReadOutputFactory> factory);
    void patch_options(boost::program_options::options_description& desc) const override;
    BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const override;
    ~OutputDispatcherFactory() override;

private:
    std::vector<std::shared_ptr<BaseReadOutputFactory>> factories_;
};

} // namespace labw::art_modern