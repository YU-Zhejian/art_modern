#include "libam_support/out/OutputDispatcher.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BamReadOutput.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/FastaReadOutput.hh"
#include "libam_support/out/FastqReadOutput.hh"
#include "libam_support/out/HeadlessBamReadOutput.hh"
#include "libam_support/out/OutParams.hh"
#include "libam_support/out/PwaReadOutput.hh"

#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>

#include <cstddef>
#include <memory>
#include <utility>

namespace labw::art_modern {

void OutputDispatcher::writeSE(const TokenRing& tokens, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    for (std::size_t i = 0; i < outputs_.size(); i++) {
        outputs_[i]->writeSE(tokens[i], pwa);
    }
}

void OutputDispatcher::writePE(const TokenRing& tokens, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    for (std::size_t i = 0; i < outputs_.size(); i++) {
        outputs_[i]->writePE(tokens[i], pwa1, pwa2);
    }
}

void OutputDispatcher::close()
{
    if (closed_) {
        return;
    }
    for (const auto& output : outputs_) {
        output->close();
    }
    closed_ = true;
}

OutputDispatcher::~OutputDispatcher() { OutputDispatcher::close(); }

void OutputDispatcher::add(std::shared_ptr<BaseReadOutput>&& output) { outputs_.emplace_back(std::move(output)); }

bool OutputDispatcher::require_alignment() const
{
    bool retv = false;

    for (const auto& output : outputs_) {
        retv |= output->require_alignment();
    }
    return retv;
}

OutputDispatcher::TokenRing OutputDispatcher::get_producer_tokens()
{
    OutputDispatcher::TokenRing retv;
    for (const auto& output : outputs_) {
        retv.emplace_back(output->get_producer_token());
    }
    return retv;
}

void OutputDispatcherFactory::patch_options(boost::program_options::options_description& desc) const
{
    for (auto const& factory : factories_) {
        factory->patch_options(desc);
    }
}
std::shared_ptr<OutputDispatcher> OutputDispatcherFactory::create(const OutParams& params) const
{
    auto output_dispatcher = std::make_shared<OutputDispatcher>();
    for (auto const& factory : factories_) {
        try {
            output_dispatcher->add(factory->create(params));
        } catch (const OutputNotSpecifiedException& e) {
            // ignored
        }
    }
    BOOST_LOG_TRIVIAL(info) << "All writers added";
    return output_dispatcher;
}

OutputDispatcherFactory::OutputDispatcherFactory()
{
    factories_.emplace_back(std::make_shared<PwaReadOutputFactory>());
    factories_.emplace_back(std::make_shared<FastaReadOutputFactory>());
    factories_.emplace_back(std::make_shared<FastqReadOutputFactory>());
    factories_.emplace_back(std::make_shared<BamReadOutputFactory>());
    factories_.emplace_back(std::make_shared<HeadlessBamReadOutputFactory>());
}

OutputDispatcherFactory::~OutputDispatcherFactory() = default;

} // namespace labw::art_modern
