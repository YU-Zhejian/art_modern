#include "bench_io_funcs.hh"

#include <mpi.h>

#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

static void mpi_job()
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    const auto processor_name_str = std::string(processor_name, name_len);

    std::ostringstream ss;
    ss << "test-" << rank;
    const std::string fpath = ss.str();
    const auto block_size = 4ULL * K_SIZE;
    const auto n_blocks = M_SIZE / size;
    const auto s_seq_write = speed_seq_write(fpath, block_size, n_blocks);
    const auto s_seq_read = speed_seq_read(fpath, block_size, n_blocks);
    const auto s_rand_read = speed_rand_read(fpath, block_size * n_blocks, block_size, n_blocks);
    { // This stuff should have been protected by an mutex.
        std::cout << "Rank " << rank << " on: " << processor_name_str << std::endl;
        std::cout << "Rank " << rank << " seq_write: " << s_seq_write / M_SIZE << " MB/s" << std::endl;
        std::cout << "Rank " << rank << " seq_read: " << s_seq_read / M_SIZE << " MB/s" << std::endl;
        std::cout << "Rank " << rank << " rand_read: " << s_rand_read / M_SIZE << " MB/s" << std::endl;
    }
    std::remove(fpath.c_str());
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    mpi_job();
    MPI_Finalize();
}
