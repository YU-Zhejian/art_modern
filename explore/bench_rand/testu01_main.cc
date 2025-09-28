// TestU01 library does not have C++ headers
extern "C" {
#include <TestU01.h>
}

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int main()
{
    swrite_Basic = FALSE; // Disable detailed TestU01 output to stdout
    std::string const temp_file = "/dev/stdin";
    char* filename = new char[temp_file.size() + 1];
    std::strncpy(filename, temp_file.c_str(), temp_file.size() + 1);
    unif01_Gen* gen = ufile_CreateReadBin(filename, 4096);
    delete[] filename;
    bbattery_Crush(gen);
    unif01_DeleteExternGenBits(gen);
    return EXIT_SUCCESS;
}
