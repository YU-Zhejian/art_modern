#include "libam_support/out/PwaReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/DumbReadOutput.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <boost/algorithm/string/join.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

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

void PwaReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    auto os = std::make_unique<std::string>(pwa.serialize());
    lfio_.push(std::move(os));
}

void PwaReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    auto os1 = std::make_unique<std::string>(pwa1.serialize());
    lfio_.push(std::move(os1));
    auto os2 = std::make_unique<std::string>(pwa2.serialize());
    lfio_.push(std::move(os2));
}

PwaReadOutput::~PwaReadOutput() { PwaReadOutput::close(); }
PwaReadOutput::PwaReadOutput(const std::string& filename, const std::vector<std::string>& args)
    : lfio_("PWA", filename, preamble(args))
{
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

void PwaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description pwa_desc("PWA Output");
    pwa_desc.add_options()("o-pwa", boost::program_options::value<std::string>(),
        "Destination of output pwa file. Unset to disable the writer.");
    desc.add(pwa_desc);
}
std::shared_ptr<BaseReadOutput> PwaReadOutputFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* /*fasta_fetch*/, const std::vector<std::string>& args) const
{
    if (vm.count("o-pwa") != 0U) {
        return std::make_shared<PwaReadOutput>(vm["o-pwa"].as<std::string>(), args);
    }
    return std::make_shared<DumbReadOutput>();
}
} // namespace labw::art_modern