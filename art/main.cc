#include <boost/log/trivial.hpp>
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include "ArtParams.hh"
#include "main_fn.hh"
#include "global_variables.hh"

using namespace std;
using namespace labw::art_modern;

int main(int argc, char* argv[])
{
#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::auto_cpu_timer t;
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif
    for (auto i = 0; i < argc; i++) {
        args.emplace_back(argv[i]);
    }
    print_banner();
    ArtParams art_params;
    art_params.parse_args(args);
    art_params.validate_args();
    art_params.print_params();
    generate_all(art_params);
    return EXIT_SUCCESS;
}
