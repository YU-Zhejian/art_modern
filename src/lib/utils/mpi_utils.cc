#include "mpi_utils.hh"

#include "art_modern_config.h" // NOLINT: for WITH_MPI
#include "art_modern_constants.hh" // NOLINT: Used in MPI functions

#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <boost/log/trivial.hpp>

#include <cstdlib>
#include <cstring> // NOLINT: Used in MPI functions
#include <string>

namespace labw::art_modern {
[[noreturn]] void exit_mpi(const int status)
{
    BOOST_LOG_TRIVIAL(info) << "EXIT";
#ifdef WITH_MPI
    int mpi_finalized_flag;
    MPI_Finalized(&mpi_finalized_flag);
    if (!mpi_finalized_flag) {
        BOOST_LOG_TRIVIAL(debug) << "Finalizing MPI...";
        MPI_Finalize();
        BOOST_LOG_TRIVIAL(debug) << "MPI finalized.";
    } else {
        BOOST_LOG_TRIVIAL(debug) << "MPI already finalized.";
    }
#endif
    std::exit(status);
}

[[noreturn]] void abort_mpi([[maybe_unused]] const int status)
{
    BOOST_LOG_TRIVIAL(info) << "ABORT";
#ifdef WITH_MPI
    BOOST_LOG_TRIVIAL(debug) << "Sending MPI_ABORT...";
    MPI_Abort(MPI_COMM_WORLD, status);
#endif
    BOOST_LOG_TRIVIAL(debug) << "Sending std::abort...";
    std::abort();
}

void bye_mpi()
{
#ifdef WITH_MPI
    BOOST_LOG_TRIVIAL(info) << "Sending BYE signal...";
    auto bufflen = static_cast<int>(std::strlen(MPI_MESSAGE_BYE));
    void* buffer = std::calloc(bufflen, sizeof(char));
    std::memcpy(buffer, MPI_MESSAGE_BYE, bufflen);
    MPI_Request request;

    MPI_Ibcast(buffer, bufflen, MPI_CHAR, MPI_MAIN_RANK, MPI_COMM_WORLD, &request);
    std::free(buffer);
#endif
}
void init_mpi([[maybe_unused]] int* argc, [[maybe_unused]] char*** argv)
{
#ifdef WITH_MPI
    MPI_Init(argc, argv);
#endif
}
std::string mpi_rank()
{
#ifdef WITH_MPI
    int mpi_finalized_flag;
    MPI_Finalized(&mpi_finalized_flag);
    if (mpi_finalized_flag) {
        return std::to_string(MPI_UNAVAILABLE_RANK);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return std::to_string(rank);
#else
    return "nompi";
#endif
}
} // namespace labw::art_modern
