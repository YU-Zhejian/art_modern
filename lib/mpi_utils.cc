#include "mpi_utils.hh"
#include "art_modern_config.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace labw {
namespace art_modern {
    void exit_mpi(const int status)
    {
        if (status == EXIT_SUCCESS) {
#ifdef WITH_MPI
            int mpi_finalized_flag;
            MPI_Finalized(&mpi_finalized_flag);
            if (!mpi_finalized_flag) {
                MPI_Finalize();
            }
#endif
            exit(status);
        } else {
            abort_mpi(status);
        }
    }

    void abort_mpi(const int status)
    {
#ifdef WITH_MPI
        MPI_Abort(MPI_COMM_WORLD, status);
#else
        std::abort();
#endif
    }
}
}
