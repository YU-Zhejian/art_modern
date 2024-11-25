#include <boost/log/trivial.hpp>

#include "DumbReadOutput.hh"
#include "FastqReadOutput.hh"
namespace labw::art_modern {

void FastqReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_ << "@" << pwa.read_name << "\n" << pwa.query << "\n+\n" << pwa.qual << "\n";
}

void FastqReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_ << "@" << pwa1.read_name << "/1\n" << pwa1.query << "\n+\n" << pwa1.qual << "\n";
    file_ << "@" << pwa2.read_name << "/2\n" << pwa2.query << "\n+\n" << pwa2.qual << "\n";
}

FastqReadOutput::~FastqReadOutput() { FastqReadOutput::close(); }
FastqReadOutput::FastqReadOutput(const std::string& filename)
    : file_(filename),
      is_closed_(false)
{
    BOOST_LOG_TRIVIAL(info) << "Writer to '" << filename << "' added.";
}

void FastqReadOutput::close() {
    if (is_closed_) {
        return;
    }
    std::scoped_lock lock(mutex_);
    file_.flush();
    file_.close();
    is_closed_ = true;
}
void FastqReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description fastq_desc("FASTQ Output");
    fastq_desc.add_options()("o-fastq", boost::program_options::value<std::string>(),
        "Destination of output FASTQ file. Unset to disable the writer.");
    desc.add(fastq_desc);
}
BaseReadOutput* FastqReadOutputFactory::create(
    const boost::program_options::variables_map& vm, const BaseFastaFetch*, const std::vector<std::string>& args) const
{
    if (vm.count("o-fastq")) {
        return new FastqReadOutput(vm["o-fastq"].as<std::string>());
    }
    return new DumbReadOutput();
}
}