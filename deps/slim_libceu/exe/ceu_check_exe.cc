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

    std::cout << ceu_interpret_cxx_std_version();
    std::cout << ceu_interpret_cxx_stdlib();
    std::cout << ceu_check_get_cxx_compiler_info();
    std::cout << ceu_check_get_ctypes_limit_info();
    std::cout << ceu_check_get_compile_time_os_info();
    std::cout << ceu_check_get_run_time_os_info();
    return EXIT_SUCCESS;
}
