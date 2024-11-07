#pragma once
#include <boost/program_options/options_description.hpp>

#include "ArtParams.hh"
#include "out/OutputDispatcher.hh"

namespace labw::art_modern {
class ArtCmdException : public std::exception { };
class ArtCmdNormalExit : public std::exception { };

/*!
 * @brief Parse command line arguments and return the parsed parameters.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return Parsed parameters.
 * @throws ArtCmdException if an error occurred.
 * @throws ArtCmdNormalExit if the user requested help or version.
 */
ArtParams parse_args(int argc, char** argv);

} // art_modern
// labw
