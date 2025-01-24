#include "libam/out/FastqReadOutput.hh"

#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseFileReadOutput.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/DumbReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

namespace {
    std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa)
    {
        auto outs = std::make_unique<std::string>();
        const std::size_t strsize = pwa.read_name.size() + (pwa.query.size() << 1U) + 6;
        outs->resize(strsize);
        outs->at(0) = 0;
        std::snprintf(
            outs->data(), strsize + 1, "@%s\n%s\n+\n%s\n", pwa.read_name.c_str(), pwa.query.c_str(), pwa.qual.c_str());
        return outs;
    }

    std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        auto outs = std::make_unique<std::string>();
        const std::size_t strsize = (pwa1.read_name.size() + (pwa1.query.size() << 1U) + 8)
            + (pwa2.read_name.size() + (pwa2.query.size() << 1U) + 8);
        outs->resize(strsize);
        outs->at(0) = 0;
        std::snprintf(outs->data(), strsize + 1, "@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", pwa1.read_name.c_str(),
            pwa1.query.c_str(), pwa1.qual.c_str(), pwa2.read_name.c_str(), pwa2.query.c_str(), pwa2.qual.c_str());
        return outs;
    }

} // namespace

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
    lfio_.push(format_fastq(pwa1, pwa2));
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
std::shared_ptr<BaseReadOutput> FastqReadOutputFactory::create(const boost::program_options::variables_map& vm,
    [[maybe_unused]] const BaseFastaFetch* /*fasta_fetch*/,
    [[maybe_unused]] const std::vector<std::string>& /*args*/) const
{
    if (vm.count("o-fastq") != 0) {
        return std::make_shared<FastqReadOutput>(vm["o-fastq"].as<std::string>());
    }
    return std::make_shared<DumbReadOutput>();
}
} // namespace labw::art_modern
