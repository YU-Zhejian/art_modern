/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/
#include "art_modern_config.h" // NOLINT: For WITH_MPI

#include "SeedAlloc.hh"
#include "libam_support/Dtypes.h"

#include "libam_support/Constants.hh"
#include "libam_support/ds/pcg_32_c.hh"
#include "libam_support/utils/mpi_utils.hh"

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <boost/log/trivial.hpp>

#include <cstdlib>
#include <ios>
#include <memory>

namespace labw::art_modern {
void SeedAlloc::run_seedalloc(const am_rand_seed_t seed)
{
    // Generate master seed
    // Algorithm: If on MPI main process, use seed from command line (if not ignored) or generate a random seed.
    // Otherwise, receive the seed from the main process.

    if (is_on_mpi_main_process_or_nompi()) {
        master_seed_ = seed;
        BOOST_LOG_TRIVIAL(debug) << "Using seed: " << std::hex << "0x" << master_seed_ << ".";
#ifdef WITH_MPI
        // Broadcast master seed to other processes
        const am_mpi_size_t size = mpi_size();
        for (am_mpi_rank_t rank = 0; rank < size; ++rank) {
            if (rank != MPI_MAIN_RANK) {
                BOOST_LOG_TRIVIAL(debug) << "Sending master seed to process " << rank << ".";
                MPI_Send(&master_seed_, 1, MPI_UINT64_T, rank, MASTER_SEED_MPI_TAG, MPI_COMM_WORLD);
            }
        }
#endif
    } else {
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
    pcg_rng_ = std::make_unique<pcg32_c>(this_process_master_seed_);
}
am_rand_seed_t SeedAlloc::nextseed() const { return pcg_rng_->operator()(); }
} // namespace labw::art_modern
