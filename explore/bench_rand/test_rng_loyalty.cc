#include "rprobs.hh"

#include <cstring>
#include <fstream>
#include <iostream>
#include <random>

namespace {

std::string to_hex(char* buff)
{
    std::string res;
    for (int i = 0; i < 4; i++) {
        char byte = buff[i];
        for (int j = 0; j < 2; j++) {
            char nibble = (byte >> ((1 - j) * 4)) & 0x0F;
            if (nibble < 10) {
                res += ('0' + nibble);
            } else {
                res += ('a' + (nibble - 10));
            }
        }
    }
    return res;
}
} // namespace
int main()
{
    std::size_t i = 0;
    char buffer[4];
    char reconst_buffer[4];
    CustomRandomDevice rng { "test.bin" };
    std::ifstream f { "test.bin", std::ios::binary };
    std::uniform_int_distribution<uint32_t> dist { std::numeric_limits<uint32_t>::min(),
        std::numeric_limits<uint32_t>::max() };
    while (f.read(buffer, sizeof(buffer))) {
        auto n = dist(rng);
        std::memcpy(reconst_buffer, &n, sizeof(n));
        if (std::strncmp(buffer, reconst_buffer, sizeof(buffer)) != 0) {
            std::cout << "Mismatch at " << i << ": " << "rng=" << to_hex(reconst_buffer) << " f=" << to_hex(buffer)
                      << "\n";
            return EXIT_FAILURE;
        }
        if (i % 1000000 == 0) {
            std::cout << "Match at " << i << ": " << "rng=" << to_hex(reconst_buffer) << " f=" << to_hex(buffer)
                      << "\n";
        }
        i++;
    }
}
