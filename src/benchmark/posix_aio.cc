// Benchmark POSIX AIO using X writes on blocks sized Y.

#include "benchmark_utils.hh"

#include <fcntl.h>
#include <unistd.h>
#include <aio.h>

#include <cerrno>
#include <chrono>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <thread>

namespace{

    constexpr std::size_t N_REPLICA = 20ULL;

void run_aio(const int fd, const size_t block_size, const size_t num_blocks) {
  std::vector<aiocb> aiocbs(num_blocks);

  for (size_t i = 0; i < num_blocks; ++i) {
    aiocbs[i].aio_fildes = fd;
    aiocbs[i].aio_buf = std::malloc(block_size);
    aiocbs[i].aio_nbytes = block_size;
  }
    for (size_t i = 0; i < num_blocks; ++i) {
        if (aio_write(&aiocbs[i]) == -1) {
        std::cerr << "Error initiating aio_write: " << std::strerror(errno) << std::endl;
        return;
        }
    }
    for (size_t i = 0; i < num_blocks; ++i) {
        while (aio_error(&aiocbs[i]) == EINPROGRESS) {
            std::this_thread::yield(); // Yield to avoid busy waiting
        }
        int ret = aio_return(&aiocbs[i]);
        if (ret == -1) {
            std::cerr << "Error in aio_return: " << std::strerror(errno) << std::endl;
        }
    }
    for (size_t i = 0; i < num_blocks; ++i) {
        std::free(aiocbs[i].aio_buf);
    }
    std::cout << std::setw(10) << block_size << " " << std::setw(10) << num_blocks << " " describe(time) << " us" << std::endl;
}

void bench(const std::string& path, const size_t block_size, const size_t num_blocks){
    int fd = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC | O_APPEND, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        std::cerr << "Error opening file: " << std::strerror(errno) << std::endl;
        return;
    }
    std::vector<std::size_t> times;
    for (auto i = 0; i < N_REPLICA; ++i){
        auto start = std::chrono::high_resolution_clock::now();
        run_aio(fd, block_size, num_blocks);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        times.emplace_back(duration);
    }
}

int main(int argc, char* argv[]) {
    std::string path = "bench.bin";
    bench(path, 256, K_SIZE * 3200ULL); // 3200 blocks of 256
    bench(path, 512, K_SIZE * 1600ULL); // 1600 blocks of 512
    bench(path, K_SIZE, K_SIZE * 800ULL); // 800 blocks of 1K
    bench(path, K_SIZE * 2ULL, K_SIZE * 400ULL); // 400 blocks of 2K
    bench(path, K_SIZE * 4ULL, K_SIZE * 200ULL); // 200 blocks of 4K
    bench(path, K_SIZE * 8ULL, K_SIZE * 100ULL); // 100 blocks of 8K
    bench(path, K_SIZE * 16ULL, K_SIZE * 50ULL); // 50 blocks of 16K
    bench(path, K_SIZE * 32ULL, K_SIZE * 25ULL); // 25 blocks of 32K
    return 0;
}

} // namespace