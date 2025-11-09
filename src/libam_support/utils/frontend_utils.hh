#pragma once

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include <cstdlib>
#include <string>

namespace labw::art_modern {
constexpr char ARG_VERSION[] = "version";
constexpr char ARG_HELP[] = "help";
void print_help(
        const boost::program_options::options_description& po_desc,
        const std::string& prefix,
        const std::string& suffix);

boost::program_options::variables_map generate_vm_while_handling_help_version(
    const boost::program_options::options_description& po_desc, int argc, char** argv,
    const std::string& prefix = "",
    const std::string& suffix = ""
    );

boost::program_options::options_description general_options();
/**
 * Extract n_threads from parallel argument.
 */
std::size_t n_threads_from_parallel(int parallel_arg);
} // namespace labw::art_modern
