#include "mpi_utils.hh"
#include "art_modern_config.h"
#ifdef WITH_MPI
#include <cstring>
#include <mpi.h>
#endif

namespace labw::art_modern {
[[noreturn]] void exit_mpi(const int status)
{
#ifdef WITH_MPI
    int mpi_finalized_flag;
    MPI_Finalized(&mpi_finalized_flag);
    if (!mpi_finalized_flag) {
        MPI_Finalize();
    }
#endif
    exit(status);
}

[[noreturn]] void abort_mpi(const int status)
{
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD, status);
#else
    std::abort();
#endif
}

void bye_mpi()
{
#ifdef WITH_MPI
    auto bufflen = static_cast<int>(std::strlen(MPI_MESSAGE_BYE));
    void* buffer = std::malloc(bufflen);
    std::memcpy(buffer, MPI_MESSAGE_BYE, bufflen);
    MPI_Bcast(buffer, bufflen, MPI_CHAR, 0, MPI_COMM_WORLD);
    std::free(buffer);
#endif
}

}
