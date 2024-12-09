#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>

#include "DumbReadOutput.hh"
#include "PwaReadOutput.hh"
namespace labw::art_modern {

void PwaReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    auto os = std::unique_ptr<std::ostringstream>();
    pwa.serialize(*os);

    lfio_.push(std::move(os));
}

void PwaReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    auto os = std::make_unique <std::ostringstream>();
    pwa1.serialize(*os);
    pwa2.serialize(*os);

    lfio_.push(std::move(os));
}

PwaReadOutput::~PwaReadOutput() { PwaReadOutput::close(); }
PwaReadOutput::PwaReadOutput(const std::string& filename, const std::vector<std::string>& args)
    : file_(filename)
    , BaseFileReadOutput(filename)
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
BaseReadOutput* PwaReadOutputFactory::create(
    const boost::program_options::variables_map& vm, const BaseFastaFetch*, const std::vector<std::string>& args) const
{
    if (vm.count("o-pwa")) {
        return new PwaReadOutput(vm["o-pwa"].as<std::string>(), args);
    }
    return new DumbReadOutput();
}
}