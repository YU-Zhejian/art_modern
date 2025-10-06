// Benchmark POSIX AIO using X writes on blocks sized Y.

#include "benchmark_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/utils/si_utils.hh"

#include <aio.h>
#include <fcntl.h>
#include <iomanip>
#include <unistd.h>

#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace labw::art_modern;

namespace {

constexpr std::size_t N_REPLICA = 20ULL;

std::size_t run_lio(const std::string& path, const size_t block_size, const size_t num_blocks)
{
    auto const fd = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC | O_APPEND | O_CLOEXEC, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        std::cerr << "Error opening file: " << std::strerror(errno) << std::endl;
        return 0;
    }
    auto* aiocbs = static_cast<aiocb **>(std::calloc(num_blocks, sizeof(aiocb*)));

    for (size_t i = 0; i < num_blocks; ++i) {
        aiocbs[i] = static_cast<aiocb*>(std::calloc(1, sizeof(aiocb)));
        aiocbs[i]->aio_fildes = fd;
        aiocbs[i]->aio_buf = std::malloc(block_size);
        aiocbs[i]->aio_nbytes = block_size;
        aiocbs[i]->aio_lio_opcode = LIO_WRITE;
    }

    auto start = std::chrono::high_resolution_clock::now();
    if (lio_listio( LIO_WAIT, aiocbs, num_blocks, nullptr)== -1) {
        std::cerr << "Error initiating lio_listio: " << std::strerror(errno) << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();


    for (size_t i = 0; i < num_blocks; ++i) {
        std::free(const_cast<void*>(aiocbs[i]->aio_buf));
        std::free(aiocbs[i]);
    }
    std::free(aiocbs);
    unlink(path.c_str());
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

std::size_t run_aio(const std::string& path, const size_t block_size, const size_t num_blocks)
{
    auto const fd = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC | O_APPEND | O_CLOEXEC, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        std::cerr << "Error opening file: " << std::strerror(errno) << std::endl;
        return 0 ;
    }
    std::vector<aiocb> aiocbs(num_blocks);

    for (size_t i = 0; i < num_blocks; ++i) {
        aiocbs[i].aio_fildes = fd;
        aiocbs[i].aio_buf = std::malloc(block_size);
        aiocbs[i].aio_nbytes = block_size;
    }
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_blocks; ++i) {
        if (aio_write(&aiocbs[i]) == -1) {
            std::cerr << "Error initiating aio_write: " << std::strerror(errno) << std::endl;
            return 0 ;
        }
    }
    for (size_t i = 0; i < num_blocks; ++i) {
        while (aio_error(&aiocbs[i]) == EINPROGRESS) {
            std::this_thread::sleep_for( std::chrono::microseconds(50) );
        }
        if (aio_return(&aiocbs[i]) == -1) {
            std::cerr << "Error in aio_return: " << std::strerror(errno) << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_blocks; ++i) {
        std::free(const_cast<void*>(aiocbs[i].aio_buf));
    }
    unlink(path.c_str());
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

void bench(const std::string& path, const size_t block_size, const size_t num_blocks)
{
    std::vector<std::size_t> times;
    for (std::size_t i = 0; i < N_REPLICA; ++i) {
        times.emplace_back(run_aio(path, block_size, num_blocks));
    }
    std::cout << "AIO: " << std::setw(10) << to_si(block_size) << "B " << std::setw(10) << num_blocks << " " << describe(times)
              << " us; mean speed=" << to_si(block_size * num_blocks / mean(times) * 1000000) << "B/s" << std::endl;

    times.clear();

    for (std::size_t i = 0; i < N_REPLICA; ++i) {
        times.emplace_back(run_lio(path, block_size, num_blocks));
    }
    std::cout << "LIO: " << std::setw(10) << to_si(block_size) << "B " << std::setw(10) << num_blocks << " " << describe(times)
              << " us; mean speed=" << to_si(block_size * num_blocks / mean(times) * 1000000) << "B/s" << std::endl;
}

} // namespace

int main()
{
    const std::string path = "out.bin";
    for (int i = 0; i < 8; i += 2) {
        bench(path, K_SIZE << i, K_SIZE >> i);
    }
    for (int i = 0; i < 8; i += 2) {
        bench(path, M_SIZE << i, K_SIZE >> i);
    }
    return 0;
}
