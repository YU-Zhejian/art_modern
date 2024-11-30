#pragma once
#include "PairwiseAlignment.hh"
#include "fasta/BaseFastaFetch.hh"
#include <boost/program_options.hpp>
#include <string>

namespace labw::art_modern {

class BaseReadOutput {
public:
    BaseReadOutput(BaseReadOutput&& other) = delete;
    BaseReadOutput(const BaseReadOutput&) = delete;
    BaseReadOutput& operator=(BaseReadOutput&&) = delete;
    BaseReadOutput& operator=(const BaseReadOutput&) = delete;

    BaseReadOutput() = default;
    virtual void writeSE(const PairwiseAlignment& pwa) = 0;
    virtual void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) = 0;
    virtual void close() = 0;
    virtual ~BaseReadOutput();
};

class BaseReadOutputFactory {
public:
    virtual const std::string name() const = 0;
    virtual void patch_options(boost::program_options::options_description& desc) const = 0;
    virtual BaseReadOutput* create(const boost::program_options::variables_map& vm, const BaseFastaFetch* fasta_fetch,
        const std::vector<std::string>& args) const
        = 0;
    virtual ~BaseReadOutputFactory();
};

} // namespace labw::art_modern
