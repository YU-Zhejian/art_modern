#include "libam/out/OutputDispatcher.hh"

#include "BaseReadOutput.hh"
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BamReadOutput.hh"
#include "libam/out/FastqReadOutput.hh"
#include "libam/out/HeadlessBamReadOutput.hh"
#include "libam/out/PwaReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

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
    for (const auto& output : outputs_) {
        output->writeSE(pwa);
    }
}

void OutputDispatcher::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    for (const auto& output : outputs_) {
        output->writePE(pwa1, pwa2);
    }
}

void OutputDispatcher::close()
{
    for (const auto& output : outputs_) {
        output->close();
    }
}

OutputDispatcher::~OutputDispatcher() { OutputDispatcher::close(); }

void OutputDispatcher::add(std::shared_ptr<BaseReadOutput>&& output) { outputs_.emplace_back(std::move(output)); }

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
    add(std::make_shared<FastqReadOutputFactory>());
    add(std::make_shared<BamReadOutputFactory>());
    add(std::make_shared<HeadlessBamReadOutputFactory>());
}

OutputDispatcherFactory::~OutputDispatcherFactory() = default;

} // namespace labw::art_modern
