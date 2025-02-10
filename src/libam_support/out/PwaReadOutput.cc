#include "libam_support/out/PwaReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"

#include <concurrentqueue.h>

#include <boost/algorithm/string/join.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {
namespace {
    std::string preamble(const std::vector<std::string>& args)
    {
        std::ostringstream oss;
        oss << "#PWA\n";
        oss << "#ARGS: " << boost::algorithm::join(args, " ") << "\n";
        return oss.str();
    }
} // namespace

void PwaReadOutput::writeSE(const moodycamel::ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    auto os = std::make_unique<std::string>(pwa.serialize());
    lfio_.push(std::move(os), token);
}

void PwaReadOutput::writePE(const moodycamel::ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    auto os1 = std::make_unique<std::string>(pwa1.serialize());
    lfio_.push(std::move(os1), token);
    auto os2 = std::make_unique<std::string>(pwa2.serialize());
    lfio_.push(std::move(os2), token);
}

PwaReadOutput::~PwaReadOutput() { PwaReadOutput::close(); }
PwaReadOutput::PwaReadOutput(const std::string& filename, const std::vector<std::string>& args, const int n_threads)
    : lfio_("PWA", filename, preamble(args))
{
    lfio_.init_queue(n_threads, 0);
    lfio_.start();
}

void PwaReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    closed_ = true;
}

bool PwaReadOutput::require_alignment() const { return true; }

    moodycamel::ProducerToken PwaReadOutput::get_producer_token() {
        return lfio_.get_producer_token();
    }

    void PwaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description pwa_desc("PWA Output");
    pwa_desc.add_options()("o-pwa", boost::program_options::value<std::string>(),
        "Destination of output pwa file. Unset to disable the writer.");
    desc.add(pwa_desc);
}
std::shared_ptr<BaseReadOutput> PwaReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-pwa") != 0U) {
        return std::make_shared<PwaReadOutput>(params.vm["o-pwa"].as<std::string>(), params.args, params.n_threads);
    }
    throw OutputNotSpecifiedException{};
}
} // namespace labw::art_modern