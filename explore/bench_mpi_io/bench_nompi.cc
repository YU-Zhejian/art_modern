#include "bench_io_funcs.hh"

#include <cstdio>
#include <iostream>
#include <string>

int main()
{
    const std::string fpath = "test";
    const auto block_size = 4ULL * K_SIZE;
    const auto n_blocks = M_SIZE;
    const auto s_seq_write = speed_seq_write(fpath, block_size, n_blocks);
    const auto s_seq_read = speed_seq_read(fpath, block_size, n_blocks);
    const auto s_rand_read = speed_rand_read(fpath, block_size * n_blocks, block_size, n_blocks);

    std::cout << "seq_write: " << s_seq_write / M_SIZE << " MB/s" << std::endl;
    std::cout << "seq_read: " << s_seq_read / M_SIZE << " MB/s" << std::endl;
    std::cout << "rand_read: " << s_rand_read / M_SIZE << " MB/s" << std::endl;

    std::remove(fpath.c_str());
}
