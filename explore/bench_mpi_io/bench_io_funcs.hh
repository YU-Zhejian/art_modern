#pragma once

#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

constexpr std::size_t K_SIZE = 1 << 10;
constexpr std::size_t M_SIZE = 1 << 20;
constexpr std::size_t G_SIZE = 1 << 20;
const std::string DEVNULL = "/dev/null";

static double speed_seq_write(const std::string& file_path, const std::size_t block_size, const std::size_t n_blocks)
{
    const std::string block = std::string(block_size, '\0');
    auto start = std::chrono::high_resolution_clock::now();
    auto ofs = std::ofstream(file_path, std::ios::out | std::ios::binary);
    for (std::size_t i = 0; i < n_blocks; i++) {
        ofs.write(block.data(), block_size);
    }
    ofs.close();
    auto end = std::chrono::high_resolution_clock::now();
    const double speed = 1.0 * (n_blocks * block_size)
        / std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * 1000;
    return speed;
}

static double speed_seq_read(const std::string& file_path, const std::size_t block_size, const std::size_t n_blocks)
{
    std::size_t i = 0;
    auto start = std::chrono::high_resolution_clock::now();
    auto ifs = std::ifstream(file_path, std::ios::in | std::ios::binary);
    std::vector<char> block;
    block.resize(block_size);
    while (i < n_blocks) {
        ifs.read(block.data(), block_size);
        i++;
    }
    ifs.close();
    auto end = std::chrono::high_resolution_clock::now();
    const double speed = 1.0 * (n_blocks * block_size)
        / std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * 1000;
    return speed;
}

static double speed_rand_read(const std::string& file_path, const std::size_t file_length, const std::size_t block_size,
    const std::size_t n_blocks)
{
    const std::size_t valid_end = file_length - block_size;
    std::mt19937 gen(std::random_device {}());
    std::uniform_int_distribution<std::size_t> dist(0, valid_end);
    std::size_t i = 0;
    auto start = std::chrono::high_resolution_clock::now();
    auto ifs = std::ifstream(file_path, std::ios::in | std::ios::binary);
    std::vector<char> block;
    block.resize(block_size);
    while (i < n_blocks) {
        ifs.seekg(dist(gen));
        ifs.read(block.data(), block_size);
        i++;
    }
    ifs.close();
    auto end = std::chrono::high_resolution_clock::now();
    const double speed = 1.0 * (n_blocks * block_size)
        / std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() * 1000;
    return speed;
}
