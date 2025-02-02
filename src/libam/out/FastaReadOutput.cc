#include "libam/out/FastaReadOutput.hh"

#include "libam/ds/PairwiseAlignment.hh"
#include "libam/out/BaseReadOutput.hh"
#include "libam/out/DumbReadOutput.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"

#include <fmt/core.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>

namespace labw::art_modern {

namespace {
    std::unique_ptr<std::string> format_fasta(const PairwiseAlignment& pwa)
    {
        return std::make_unique<std::string>(fmt::format(">{}\n{}\n", pwa.read_name, pwa.query));
    }

    std::unique_ptr<std::string> format_fasta(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
    {
        return std::make_unique<std::string>(
            fmt::format(">{}\n{}\n>{}\n{}\n", pwa1.read_name, pwa1.query, pwa2.read_name, pwa2.query));
    }

} // namespace

void FastaReadOutput::writeSE(const PairwiseAlignment& pwa)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fasta(pwa));
}

void FastaReadOutput::writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2)
{
    if (closed_) {
        return;
    }
    lfio_.push(format_fasta(pwa1, pwa2));
}

FastaReadOutput::~FastaReadOutput() { FastaReadOutput::close(); }
FastaReadOutput::FastaReadOutput(const std::string& filename)
    : lfio_("FASTA", filename)
{
    lfio_.start();
}

void FastaReadOutput::close()
{
    if (closed_) {
        return;
    }
    lfio_.stop();
    closed_ = true;
}

bool FastaReadOutput::require_alignment() const { return false; }

void FastaReadOutputFactory::patch_options(boost::program_options::options_description& desc) const
{
    boost::program_options::options_description fasta_desc("FASTA Output");
    fasta_desc.add_options()("o-fasta", boost::program_options::value<std::string>(),
        "Destination of output FASTA file. Unset to disable the writer.");
    desc.add(fasta_desc);
}
std::shared_ptr<BaseReadOutput> FastaReadOutputFactory::create(const boost::program_options::variables_map& vm,
    [[maybe_unused]] const BaseFastaFetch* /*fasta_fetch*/,
    [[maybe_unused]] const std::vector<std::string>& /*args*/) const
{
    if (vm.count("o-fasta") != 0) {
        return std::make_shared<FastaReadOutput>(vm["o-fasta"].as<std::string>());
    }
    return std::make_shared<DumbReadOutput>();
}
} // namespace labw::art_modern
