/**
 * @brief Argument parser
 */

#pragma once
#include "art/lib/ArtIOParams.hh"
#include "art/lib/ArtParams.hh"

#include <tuple>

namespace labw::art_modern {

/**
 * @brief Parse command line arguments and return the parsed parameters.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return Parsed parameters.
 */
std::tuple<ArtParams, ArtIOParams> parse_args(int argc, char** argv);

} // namespace labw::art_modern
