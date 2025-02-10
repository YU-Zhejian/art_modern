#pragma once

#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <boost/program_options/variables_map.hpp>

#include <memory>
#include <string>
#include <vector>

namespace labw::art_modern{

    struct OutParams{
        const int n_threads;
        const boost::program_options::variables_map vm;
        const std::vector<std::string> args;
        const std::shared_ptr<BaseFastaFetch> fasta_fetch;
    };
} // namespace labw::art_modern