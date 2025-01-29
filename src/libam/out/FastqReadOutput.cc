#include "libam/out/FastqReadOutput.hh"

#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseFileReadOutput.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/DumbReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

#include <fmt/core.h>

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
        return std::make_unique<std::string>(fmt::format("@{}\n{}\n+\n{}\n", pwa.read_name, pwa.query, pwa.qual));
    }

    std::unique_ptr<std::string> format_fastq(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        return std::make_unique<std::string>(fmt::format("@{}\n{}\n+\n{}\n@{}\n{}\n+\n{}\n", pwa1.read_name, pwa1.query,
            pwa1.qual, pwa2.read_name, pwa2.query, pwa2.qual));
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
    , lfio_("FASTQ", file_)
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

bool FastqReadOutput::require_alignment() const { return false; }

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
