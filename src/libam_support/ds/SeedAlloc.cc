#include "art_modern_config.h" // NOLINT: For WITH_MPI

#include "SeedAlloc.hh"
#include "libam_support/Dtypes.h"

#include "libam_support/Constants.hh"
#include "libam_support/ds/pcg_32_c.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/rand_utils.hh"

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <boost/log/trivial.hpp>

#include <cstdlib>
#include <vector>

namespace labw::art_modern
{
    void SeedAlloc::run_seedalloc(
        const bool ignore_seed_from_cmdline,
        const am_rand_seed_t seed_from_cmdline,
        const std::size_t n_threads
    )
        {
            // Generate master seed
            // Algorithm: If on MPI main process, use seed from command line (if not ignored) or generate a random seed.
            // Otherwise, receive the seed from the main process.

            if (is_on_mpi_main_process_or_nompi())
            {
                if (!ignore_seed_from_cmdline)
                {
                    master_seed_ = seed_from_cmdline;
                    BOOST_LOG_TRIVIAL(debug) << "Using seed from command line: " <<std::hex << "0x" <<  master_seed_ << ".";
                }
                else
                {
                    master_seed_ = rand_seed();
                    BOOST_LOG_TRIVIAL(debug) << "Generated random master seed: " <<std::hex << "0x" <<  master_seed_ << ".";
                }
#ifdef WITH_MPI
                // Broadcast master seed to other processes
                const am_mpi_size_t size = mpi_size();
                for (am_mpi_rank_t rank = 0; rank < size; ++rank)
                {
                    if (rank != MPI_MAIN_RANK)
                    {
                        BOOST_LOG_TRIVIAL(debug) << "Sending master seed to process " << rank << ".";
                        MPI_Send(&master_seed_, 1, MPI_UINT64_T, rank, MASTER_SEED_MPI_TAG, MPI_COMM_WORLD);
                    }
                }
#endif
            } else
            {
#ifdef WITH_MPI
                MPI_Status status;
                MPI_Recv(&master_seed_, 1, MPI_UINT64_T, MPI_MAIN_RANK, MASTER_SEED_MPI_TAG, MPI_COMM_WORLD, &status);
                BOOST_LOG_TRIVIAL(debug) << "Received master seed from main process: " << master_seed_ << ".";
#endif
            }
            // Use master seed shifted with MPI rank to seed
#ifdef WITH_MPI
            this_process_master_seed_ = master_seed_ + static_cast<am_rand_seed_t>(mpi_rank());
#else
            this_process_master_seed_ = master_seed_;
#endif
            // Use PCG to generate seeds
            pcg32_c pcg_rng(this_process_master_seed_);
            for (std::size_t i = 0; i < n_threads; ++i)
            {
                const auto allocated_seed = pcg_rng();
                allocated_seeds_.emplace_back(allocated_seed);
            }
        }
    am_rand_seed_t SeedAlloc::seed(const std::size_t thread_idx) const
    {
        return allocated_seeds_.at(thread_idx);
    }
} // namespace labw::art_modern
