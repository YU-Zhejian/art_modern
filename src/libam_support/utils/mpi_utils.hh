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
#include <cstdlib>
#include <string>

namespace labw::art_modern {
constexpr  char MPI_MESSAGE_BYE[] = "BYE\0";

void init_mpi(int* argc, char*** argv);

/**
 * Call `MPI_Finalize` (if MPI is not finalized) with the given status.
 *
 * Every process should call this method.
 *
 * @param status Exit status.
 */
void exit_mpi(int status);

/**
 * Call `MPI_Abort` or `std::abort` with the given status.
 * Will terminate all processes.
 *
 * @param status Exit status.
 */
[[noreturn]] void abort_mpi(int status = EXIT_FAILURE);
/**
 * Broadcast the message "BYE" to all processes.
 */
void bye_mpi();

/**
 * Get the current MPI rank in string.
 *
 * @return The current MPI rank in string.
 * MPI_UNAVAILABLE_RANK if MPI had stopped.
 * "nompi" if MPI is not available.
 */
std::string mpi_rank();

} // namespace labw::art_modern
