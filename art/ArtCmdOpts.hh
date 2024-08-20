#include <boost/program_options.hpp>
#include "ArtParams.hh"

namespace labw
{
namespace art_modern
{

boost::program_options::options_description option_parser();

class ArtCmdOpts{
    const boost::program_options::options_description po_desc = option_parser();
    ArtParams parse_args(int argc, char** argv);
};



} // art_modern
} // labw
