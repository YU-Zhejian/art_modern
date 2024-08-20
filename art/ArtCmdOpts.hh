#include <boost/program_options/options_description.hpp>

#include "ArtParams.hh"
#include "out/OutputDispatcher.hh"

namespace labw {
namespace art_modern {

    boost::program_options::options_description option_parser();
    OutputDispatcherFactory get_output_dispatcher_factory();

    class ArtCmdOpts {
    public:
        ArtParams parse_args(int argc, char** argv) const;

    private:
        const boost::program_options::options_description po_desc_ = option_parser();
        const OutputDispatcherFactory out_dispatcher_factory_ = get_output_dispatcher_factory();
    };

} // art_modern
} // labw
