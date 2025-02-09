#include "libam_support/out/OutputDispatcher.hh"

#include "BaseReadOutput.hh"
#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BamReadOutput.hh"
#include "libam_support/out/FastaReadOutput.hh"
#include "libam_support/out/FastqReadOutput.hh"
#include "libam_support/out/HeadlessBamReadOutput.hh"
#include "libam_support/out/PwaReadOutput.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <boost/log/trivial.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

void OutputDispatcher::writeSE(const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    for (const auto& output : outputs_) {
        output->writeSE(pwa);
    }
}

void OutputDispatcher::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    for (const auto& output : outputs_) {
        output->writePE(pwa1, pwa2);
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

void OutputDispatcherFactory::patch_options(boost::program_options::options_description& desc) const
{
    for (auto const& factory : factories_) {
        factory->patch_options(desc);
    }
}
std::shared_ptr<BaseReadOutput> OutputDispatcherFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const
{
    auto output_dispatcher = std::make_shared<OutputDispatcher>();
    for (auto const& factory : factories_) {
        output_dispatcher->add(factory->create(vm, fasta_fetch, args));
    }
    BOOST_LOG_TRIVIAL(info) << "All writers added";
    return output_dispatcher;
}
void OutputDispatcherFactory::add(std::shared_ptr<BaseReadOutputFactory> factory)
{
    factories_.emplace_back(std::move(factory));
}

OutputDispatcherFactory::OutputDispatcherFactory()
{
    add(std::make_shared<PwaReadOutputFactory>());
    add(std::make_shared<FastaReadOutputFactory>());
    add(std::make_shared<FastqReadOutputFactory>());
    add(std::make_shared<BamReadOutputFactory>());
    add(std::make_shared<HeadlessBamReadOutputFactory>());
}

OutputDispatcherFactory::~OutputDispatcherFactory() = default;

} // namespace labw::art_modern
