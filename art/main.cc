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

#include "ArtCmdOpts.hh"
#include "main_fn.hh"

#include "utils/dump_utils.hh"
#include "utils/log_utils.hh"
#include "utils/mpi_utils.hh"
#include "utils/protobuf_utils.hh"

using namespace labw::art_modern;

#ifdef WITH_MPI
int main_mpi_child()
{
    char buffer[100];
    std::string received_message;
    MPI_Request request;
    MPI_Ibcast(buffer, 100, MPI_CHAR, MPI_MAIN_RANK, MPI_COMM_WORLD, &request);

    while (true) {
        BOOST_LOG_TRIVIAL(info) << "Receiving signal...";
        MPI_Status status;
        int flag;
        MPI_Test(&request, &flag, &status);
        if (flag) {
            if (std::strncmp(buffer, MPI_MESSAGE_BYE, std::strlen(MPI_MESSAGE_BYE)) == 0) {
                BOOST_LOG_TRIVIAL(info) << "Received BYE signal.";
                return EXIT_SUCCESS;
            } else {
                BOOST_LOG_TRIVIAL(info) << "Received wrong signal.";
                sleep(1);
            }
        } else {
            BOOST_LOG_TRIVIAL(info) << "Received no signal.";
            sleep(1);
        }
    }
}

#endif

void handle_mpi_child()
{
#ifdef WITH_MPI
    int mpi_comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_comm_rank);
    if (mpi_comm_rank == MPI_MAIN_RANK) {
        print_banner();
        int mpi_comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
        BOOST_LOG_TRIVIAL(info) << "MPI detected with " << mpi_comm_size
                                << " MPI-parallelized processes running in total.";
    } else {
        BOOST_LOG_TRIVIAL(info) << "MPI detected. This process have rank " << mpi_comm_rank << ".";
        exit_mpi(main_mpi_child());
    }
#else
    BOOST_LOG_TRIVIAL(warning) << "MPI not found! Cross-node parallelism disabled.";
#endif
}

int main(int argc, char* argv[])
{
    init_mpi(&argc, &argv);
    // 1st round initialization of a working console logger
    init_logger();
    init_file_logger("log.d");
    handle_mpi_child();
    handle_dumps();
    validate_protobuf_version();

    try {
        auto parse_args_1st = parse_args(argc, argv);
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
    // FIXME: Come up a way to avoid this 2nd parsing of params.
    auto art_params = parse_args(argc, argv);
    generate_all(art_params);
#ifdef WITH_BOOST_TIMER
    t.stop();
    BOOST_LOG_TRIVIAL(info) << "Time spent: " << t.format(3, "%w");
#endif
    bye_mpi();
    exit_mpi(EXIT_SUCCESS);
}
