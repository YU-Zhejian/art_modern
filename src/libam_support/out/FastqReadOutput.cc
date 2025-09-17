#include "libam_support/out/FastqReadOutput.hh"

#include "libam_support/ds/PairwiseAlignment.hh"
#include "libam_support/lockfree/ProducerToken.hh"
#include "libam_support/out/BaseReadOutput.hh"
#include "libam_support/out/OutParams.hh"

#include <fmt/format.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

#include <cstdio>
#include <memory>
#include <string>

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

void FastqReadOutput::writeSE(const ProducerToken& token, const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fastq(pwa), token);
}

void FastqReadOutput::writePE(const ProducerToken& token, const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fastq(pwa1, pwa2), token);
}

FastqReadOutput::~FastqReadOutput() { FastqReadOutput::close(); }
FastqReadOutput::FastqReadOutput(const std::string& filename, const int n_threads)
    : lfio_("FASTQ", filename)
{
    lfio_.init_queue(n_threads, 0);
    lfio_.start();
}

void FastqReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    closed_ = true;
}

bool FastqReadOutput::require_alignment() const { return false; }

ProducerToken FastqReadOutput::get_producer_token() { return lfio_.get_producer_token(); }

void FastqReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description fastq_desc("FASTQ Output");
    fastq_desc.add_options()("o-fastq", boost::program_options::value<std::string>(),
        "Destination of output FASTQ file. Unset to disable the writer.");
    desc.add(fastq_desc);
}
std::shared_ptr<BaseReadOutput> FastqReadOutputFactory::create(const OutParams& params) const
{
    if (params.vm.count("o-fastq") != 0) {
        return std::make_shared<FastqReadOutput>(params.vm["o-fastq"].as<std::string>(), params.n_threads);
    }
    throw OutputNotSpecifiedException {};
}

} // namespace labw::art_modern
