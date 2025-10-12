/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

/**
 * @brief Utility functions for MPI.
 *
 * Here defines some utility functions that would work under MPI.
 *
 */

#pragma once

#include <cstddef>
#include <cstdlib>
#include <string>

namespace labw::art_modern {

bool is_mpi_finalized();

void init_mpi(int* argc, char*** argv);

/**
 * @throw std::runtime_error if the MPI is finalized (in case init_mpi is not called).
 *
 * @return 1 if MPI is not available. Otherwise, return the size of MPI_COMM_WORLD.
 */
std::size_t mpi_size();

/**
 * Call `MPI_Finalize` (if MPI is not finalized).
 *
 * Every process should call this method.
 */
void exit_mpi();

/**
 * Call `MPI_Abort` or `std::abort` with the given status.
 * Will terminate all processes.
 *
 * @param status Exit status.
 */

[[noreturn]] void abort_mpi(int status = EXIT_FAILURE);

/**
 * Get the current MPI rank in string.
 *
 * @return The current MPI rank in string.
 * MPI_UNAVAILABLE_RANK if MPI had stopped.
 * "nompi" if MPI is not available.
 */
std::string mpi_rank();

} // namespace labw::art_modern
