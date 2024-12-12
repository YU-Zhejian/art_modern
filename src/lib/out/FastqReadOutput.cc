#include "FastqReadOutput.hh"
#include "DumbReadOutput.hh"

namespace labw::art_modern {

std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa)
{
    auto outs = std::make_unique<std::string>();
    const std::size_t strsize = pwa.read_name.size() + (pwa.query.size() << 1) + 6;
    outs->resize(strsize);
    outs->at(0) = 0;
    std::snprintf(
        outs->data(), strsize + 1, "@%s\n%s\n+\n%s\n", pwa.read_name.c_str(), pwa.query.c_str(), pwa.qual.c_str());
    return outs;
}

std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa, const bool is_read1)
{
    auto outs = std::make_unique<std::string>();
    const std::size_t strsize = pwa.read_name.size() + (pwa.query.size() << 1) + 8;
    outs->resize(strsize);
    outs->at(0) = 0;
    std::snprintf(outs->data(), strsize + 1, "@%s/%d\n%s\n+\n%s\n", pwa.read_name.c_str(), is_read1 ? 1 : 2,
        pwa.query.c_str(), pwa.qual.c_str());
    return outs;
}

void FastqReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (is_closed_) {
        return;
    }
    lfio_.push(format_fastq(pwa));
}

void FastqReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (is_closed_) {
        return;
    }
    lfio_.push(format_fastq(pwa1, true));
    lfio_.push(format_fastq(pwa2, false));
}

FastqReadOutput::~FastqReadOutput() { FastqReadOutput::close(); }
FastqReadOutput::FastqReadOutput(const std::string& filename)
    : BaseFileReadOutput(filename)
    , file_(filename)
    , lfio_(file_)
{
    lfio_.start();
}

void FastqReadOutput::close()
{
    if (is_closed_) {
        return;
    }
    lfio_.stop();
    file_.flush();
    file_.close();
    BaseFileReadOutput::close();
}
void FastqReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description fastq_desc("FASTQ Output");
    fastq_desc.add_options()("o-fastq", boost::program_options::value<std::string>(),
        "Destination of output FASTQ file. Unset to disable the writer.");
    desc.add(fastq_desc);
}
BaseReadOutput* FastqReadOutputFactory::create(const boost::program_options::variables_map& vm, const BaseFastaFetch*,
    [[maybe_unused]] const std::vector<std::string>& args) const
{
    if (vm.count("o-fastq")) {
        return new FastqReadOutput(vm["o-fastq"].as<std::string>());
    }
    return new DumbReadOutput();
}
}