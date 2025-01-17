/*!
 * @brief Utility functions for MPI.
 *
 */

#pragma once
#include <cstdlib>
#include <string>

namespace labw::art_modern {
const char MPI_MESSAGE_BYE[] = "BYE\0";

void init_mpi(int* argc, char*** argv);

/*!
 * Call `MPI_Finalize` (if MPI is not finalized) or `std::exit` with the given status.
 *
 * Every process should call this method.
 *
 * @param status Exit status.
 */
[[noreturn]] void exit_mpi(int status);

/*!
 * Call `MPI_Abort` or `std::abort` with the given status.
 * Will terminate all processes.
 *
 * @param status Exit status.
 */
[[noreturn]] void abort_mpi(int status = EXIT_FAILURE);
/*!
 * Broadcast the message "BYE" to all processes.
 */
void bye_mpi();

std::string mpi_rank();

} // namespace labw::art_modern
