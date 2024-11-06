#include "art_modern_config.h"
#include <boost/log/trivial.hpp>
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "ArtCmdOpts.hh"
#include "main_fn.hh"

int main(int argc, char* argv[])
{
#ifdef WITH_MPI
    MPI_Init(&argc, &argv);
    int mpi_comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_comm_rank);
    if (mpi_comm_rank == 0) {
        int mpi_comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
        BOOST_LOG_TRIVIAL(info) << "MPI detected with " << mpi_comm_size
                                << " MPI-parallelized processes running in total.";
    } else {
        BOOST_LOG_TRIVIAL(info) << "MPI detected. This process have rank " << mpi_comm_rank << ".";
        BOOST_LOG_TRIVIAL(info) << "MPI-based parallelization not implemented; terminating the process.";
        MPI_Finalize();
        return EXIT_SUCCESS;
    }
#else
    BOOST_LOG_TRIVIAL(warning) << "MPI not found! Cross-node parallelism disabled.";
#endif
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
#ifdef WITH_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
