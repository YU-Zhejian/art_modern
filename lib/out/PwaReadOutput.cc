#include <boost/log/trivial.hpp>

#include "DumbReadOutput.hh"
#include "PwaReadOutput.hh"
namespace labw::art_modern {

void PwaReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_ << pwa.serialize();
}

void PwaReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_ << pwa1.serialize() << pwa2.serialize();
}

PwaReadOutput::~PwaReadOutput() { PwaReadOutput::close(); }
PwaReadOutput::PwaReadOutput(const std::string& filename)
    : file_(filename)
    , filename(filename)
    , is_closed_(false)
{
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
}

void PwaReadOutput::close()
{
    if (is_closed_) {
        return;
    }
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' closed.";
    file_.flush();
    file_.close();
}
void PwaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description pwa_desc("PWA Output");
    pwa_desc.add_options()("o-pwa", boost::program_options::value<std::string>(),
        "Destination of output pwa file. Unset to disable the writer.");
    desc.add(pwa_desc);
}
BaseReadOutput* PwaReadOutputFactory::create(
    const boost::program_options::variables_map& vm, const BaseFastaFetch*, const std::vector<std::string>& args) const
{
    if (vm.count("o-pwa")) {
        return new PwaReadOutput(vm["o-pwa"].as<std::string>());
    }
    return new DumbReadOutput();
}
}