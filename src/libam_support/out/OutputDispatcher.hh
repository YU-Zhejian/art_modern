#pragma once

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>

#include <concurrentqueue.h>

#include <atomic>
#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

class OutputDispatcher {

public:
    DELETE_MOVE(OutputDispatcher)
    DELETE_COPY(OutputDispatcher)
    using TokenRing =std::vector<moodycamel::ProducerToken>;

    [[nodiscard]] bool require_alignment() const ;

    OutputDispatcher() = default;
    ~OutputDispatcher();

    void add(std::shared_ptr<BaseReadOutput>&& output);
    void writeSE(const TokenRing& tokens, const PairwiseAlignment& pwa) ;
    void writePE(const TokenRing& tokens, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) ;
    void close();
    TokenRing get_producer_tokens();

private:
    std::vector<std::shared_ptr<BaseReadOutput>> outputs_;
    std::atomic<bool> closed_ = false;
};

class OutputDispatcherFactory {
public:
    DELETE_MOVE(OutputDispatcherFactory)
    DELETE_COPY(OutputDispatcherFactory)

    OutputDispatcherFactory();

    [[nodiscard]] const std::string name() const { return "OD"; }
    void patch_options(boost::program_options::options_description& desc) const;
    [[nodiscard]] std::shared_ptr<OutputDispatcher> create(const OutParams& params) const;
    ~OutputDispatcherFactory();

private:
    std::vector<std::shared_ptr<BaseReadOutputFactory>> factories_;
};

} // namespace labw::art_modern