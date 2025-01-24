#pragma once
#include "libam/ds/PairwiseAlignment.hh"
#include "libam/ref/fetch/BaseFastaFetch.hh"
#include "libam/utils/class_macros_utils.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern {

class BaseReadOutput {
public:
    DELETE_MOVE(BaseReadOutput)
    DELETE_COPY(BaseReadOutput)

    BaseReadOutput() = default;
    virtual void writeSE(const PairwiseAlignment& pwa) = 0;
    virtual void writePE(const PairwiseAlignment& pwa1, const PairwiseAlignment& pwa2) = 0;
    virtual void close() = 0;
    virtual ~BaseReadOutput() = default;
};

class BaseReadOutputFactory {
public:
    DELETE_MOVE(BaseReadOutputFactory)
    DELETE_COPY(BaseReadOutputFactory)
    BaseReadOutputFactory() = default;

    [[nodiscard]] virtual const std::string name() const = 0;
    virtual void patch_options(boost::program_options::options_description& desc) const = 0;
    virtual std::shared_ptr<BaseReadOutput> create(const boost::program_options::variables_map& vm,
        const BaseFastaFetch* fasta_fetch, const std::vector<std::string>& args) const
        = 0;
    virtual ~BaseReadOutputFactory() = default;
};

} // namespace labw::art_modern
