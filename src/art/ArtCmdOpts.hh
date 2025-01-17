#pragma once
#include "art/ArtParams.hh"

namespace labw::art_modern {

/*!
 * @brief Parse command line arguments and return the parsed parameters.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return Parsed parameters.
 */
ArtParams parse_args(int argc, char** argv);

} // namespace labw::art_modern
