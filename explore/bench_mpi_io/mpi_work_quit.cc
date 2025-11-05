#include <chrono>
#include <cstring>
#include <iostream>
#include <mpi.h>
#include <thread>

enum Signal { WORK = 1, QUIT = 2 };

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Rank 0: decide and signal
    if (rank == 0) {
        bool do_work = false;
        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "--work") == 0) {
                do_work = true;
                break;
            }
        }
        Signal sig = do_work ? WORK : QUIT;
        for (int dest = 1; dest < size; ++dest) {
            MPI_Send(&sig, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        }
        std::cout << "[Rank 0] Sent signal: " << (sig == WORK ? "WORK" : "QUIT") << std::endl;
    } else {
        Signal sig;
        MPI_Recv(&sig, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (sig == WORK) {
            std::cout << "[Rank " << rank << "] Received WORK, sleeping 20s..." << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(20));
        } else {
            std::cout << "[Rank " << rank << "] Received QUIT, exiting." << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
