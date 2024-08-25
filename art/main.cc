#include <boost/log/trivial.hpp>
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include "ArtCmdOpts.hh"
#include "main_fn.hh"

int main(int argc, char* argv[])
{
    labw::art_modern::init_logger();
    labw::art_modern::handle_dumps();
#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::auto_cpu_timer t;
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif
    labw::art_modern::print_banner();
    labw::art_modern::ArtCmdOpts art_cmd_opts;
    auto art_params = art_cmd_opts.parse_args(argc, argv);
    BOOST_LOG_TRIVIAL(info) << "Argument parsing finished. Start generating...";
    generate_all(art_params);
    return EXIT_SUCCESS;
}
