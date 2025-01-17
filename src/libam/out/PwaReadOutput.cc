#include "libam/out/PwaReadOutput.hh"

#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseFileReadOutput.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/DumbReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

#include <boost/algorithm/string/join.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace labw::art_modern {

void PwaReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    auto os = std::make_unique<std::string>(pwa.serialize());
    lfio_.push(std::move(os));
}

void PwaReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    auto os1 = std::make_unique<std::string>(pwa1.serialize());
    lfio_.push(std::move(os1));
    auto os2 = std::make_unique<std::string>(pwa2.serialize());
    lfio_.push(std::move(os2));
}

PwaReadOutput::~PwaReadOutput() { PwaReadOutput::close(); }
PwaReadOutput::PwaReadOutput(const std::string& filename, const std::vector<std::string>& args)
    : BaseFileReadOutput(filename)
    , file_(filename)
    , lfio_(file_)
{
    file_ << "#PWA\n";
    file_ << "#ARGS: " << boost::algorithm::join(args, " ") << "\n";
    file_.flush();

    lfio_.start();
}

void PwaReadOutput::close()
{
    if (is_closed_) {
        return;
    }
    lfio_.stop();
    file_.flush();
    BaseFileReadOutput::close();
}
void PwaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description pwa_desc("PWA Output");
    pwa_desc.add_options()("o-pwa", boost::program_options::value<std::string>(),
        "Destination of output pwa file. Unset to disable the writer.");
    desc.add(pwa_desc);
}
BaseReadOutput* PwaReadOutputFactory::create(const boost::program_options::variables_map& vm,
    const BaseFastaFetch* /*fasta_fetch*/, const std::vector<std::string>& args) const
{
    if (vm.count("o-pwa") != 0U) {
        return new PwaReadOutput(vm["o-pwa"].as<std::string>(), args);
    }
    return new DumbReadOutput();
}
} // namespace labw::art_modern