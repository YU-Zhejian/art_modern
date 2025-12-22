#include "art_modern_config.h" // NOLINT: For WITH_MPI
#include "libam_support/Dtypes.h"

#include "libam_support/Constants.hh"
#include "libam_support/ds/SeedAlloc.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/seq_utils.hh"

#include <boost/log/trivial.hpp>
#include <thread>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <vector>

using namespace labw::art_modern;

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);
    init_logger();
    init_file_logger("mpi_allocseeds", true);
    std::size_t const n_threads = std::thread::hardware_concurrency();

    SeedAlloc seed_alloc;
    seed_alloc.run_seedalloc(0);
    std::vector<am_rand_seed_t> seeds;
    for (std::size_t i = 0; i < n_threads; ++i) {
        seeds.emplace_back(seed_alloc.nextseed());
    }
    // Copy them back to rank 0 process
#ifdef WITH_MPI
    if (is_on_mpi_main_process_or_nompi()) {
        const auto size = mpi_size();
        BOOST_LOG_TRIVIAL(info) << "Collecting allocated seeds from other processes.";
        BOOST_LOG_TRIVIAL(info) << "Received allocated seeds from rank " << MPI_MAIN_RANK << ": " << vec2str(seeds);
        for (am_mpi_rank_t rank = MPI_MAIN_RANK + 1; rank < size; ++rank) {
            const auto n_seeds = static_cast<int>(n_threads);
            std::vector<am_rand_seed_t> recv_seeds(n_threads);
            MPI_Recv(recv_seeds.data(), n_seeds, MPI_UINT64_T, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            BOOST_LOG_TRIVIAL(info) << "Received allocated seeds from rank " << rank << ": " << vec2str(recv_seeds);
        }

    } else {
        const auto rank = mpi_rank();
        const auto n_seeds = static_cast<int>(seeds.size());
        BOOST_LOG_TRIVIAL(info) << "Sending allocated seeds to main process from rank " << rank << ".";
        MPI_Send(seeds.data(), n_seeds, MPI_UINT64_T, MPI_MAIN_RANK, 0, MPI_COMM_WORLD);
    }
#else
    BOOST_LOG_TRIVIAL(info) << "Allocated seeds:";
#endif
    exit_mpi();
    return EXIT_SUCCESS;
}
