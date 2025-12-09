#include "ceu_check/ceu_check_arch.h"
#include "ceu_check/ceu_check_arch_macros.h"
#include "ceu_check/ceu_check_c_compiler.h"
#include "ceu_check/ceu_check_c_std.h"
#include "ceu_check/ceu_check_c_stdlib.h"

#include "ceu_check/ceu_check_ctypes_limit.hh"
#include "ceu_check/ceu_check_cxx_compiler.hh"
#include "ceu_check/ceu_check_cxx_std.hh"
#include "ceu_check/ceu_check_cxx_stdlib.hh"
#include "ceu_check/ceu_check_os.hh"

#include <cstdlib>
#include <iostream>

int main()
{
    {
        std::cout << "Architecture: " << CEU_ARCH_NAME << std::endl;
        std::cout << "\tEndianness: ";
#if defined(CEU_COMPILE_TIME_IS_LITTLE_ENDIAN)
        std::cout << "Little Endian (Compile-time)";
#elif defined(CEU_COMPILE_TIME_IS_BIG_ENDIAN)
        std::cout << "Big Endian (Compile-time)";
#else
        if (ceu_is_little_endian()) {
            std::cout << "Little Endian (Run-time)";
        } else if (ceu_is_big_endian()) {
            std::cout << "Big Endian (Run-time)";
        } else {
            std::cout << "Unknown (Run-time)";
        }
#endif
        std::cout << std::endl;
        std::cout << "\tBitness: ";
#if defined(CEU_ARCHITECTURE_16_BIT)
        std::cout << "16-bit (Compile-time Macro)";
#elif defined(CEU_ARCHITECTURE_32_BIT)
        std::cout << "32-bit (Compile-time Macro)";
#elif defined(CEU_ARCHITECTURE_64_BIT)
        std::cout << "64-bit (Compile-time Macro)";
#else
        if (ceu_is_system_16_bit()) {
            std::cout << "16-bit (Compile-time Sizeof)";
        } else if (ceu_is_system_32_bit()) {
            std::cout << "32-bit (Compile-time Sizeof)";
        } else if (ceu_is_system_64_bit()) {
            std::cout << "64-bit (Compile-time Sizeof)";
        } else {
            std::cout << "Unknown (Compile-time Sizeof)";
        }
#endif
        std::cout << std::endl;
    }
    std::cout << std::endl;
    {
        char* info = ceu_interpret_c_std_version();
        if (info == nullptr) {
            std::cerr << "Failed to get C standard information." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << info;
        std::free(info);
    }
    {
        char* info = ceu_interpret_c_stdlib_version();
        if (info == nullptr) {
            std::cerr << "Failed to get C standard library version information." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << info;
        std::free(info);
    }
    {
        char* info = ceu_check_get_c_compiler_info();
        if (info == nullptr) {
            std::cerr << "Failed to get C compiler information." << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << info;
        std::free(info);
    }
    std::cout << std::endl;

    std::cout << ceu_interpret_cxx_std_version();
    std::cout << ceu_interpret_cxx_stdlib();
    std::cout << ceu_check_get_cxx_compiler_info();
    std::cout << std::endl;
    std::cout << ceu_check_get_ctypes_limit_info();
    std::cout << std::endl;
    std::cout << ceu_check_get_compile_time_os_info();
    std::cout << ceu_check_get_run_time_os_info();
    return EXIT_SUCCESS;
}
