#include "art_modern_config.h"
#include <boost/log/trivial.hpp>

// Boost timer
#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

// MPI
#ifdef WITH_MPI
#include <mpi.h>
#endif

// Protobuf
#ifdef WITH_PROTOBUF
#include <google/protobuf/stubs/common.h>
#endif

#include "ArtCmdOpts.hh"
#include "main_fn.hh"

#include "utils/mpi_utils.hh"

using namespace labw::art_modern;

int main(int argc, char* argv[])
{
#ifdef WITH_PROTOBUF
    GOOGLE_PROTOBUF_VERIFY_VERSION;
#endif

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
    init_logger();
    handle_dumps();
    print_banner();
    try {
        auto art_params = parse_args(argc, argv);
    } catch (const ArtCmdNormalExit&) {
        bye_mpi();
        exit_mpi(EXIT_SUCCESS);
    } catch (const ArtCmdException&) {
        abort_mpi(EXIT_FAILURE);
    }
    BOOST_LOG_TRIVIAL(info) << "Argument parsing finished. Start generating...";

#ifdef WITH_BOOST_TIMER
    BOOST_LOG_TRIVIAL(info) << "Boost::timer started.";
    boost::timer::cpu_timer t;
    t.start();
#else
    BOOST_LOG_TRIVIAL(warning) << "Boost::timer not found! Resource consumption statistics disabled.";
#endif
    // generate_all(art_params); // FIXME
#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time spent: " << t.format(3, "%w");
#endif
    bye_mpi();
    exit_mpi(EXIT_SUCCESS);
}
