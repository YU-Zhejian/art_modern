#include "ceu_check/ceu_check_cxx_std.hh"
#include "ceu_check/ceu_check_cc.hh"
#include "ceu_check/ceu_check_ctypes_limit.hh"
#include "ceu_check/ceu_check_os.hh"
#include "ceu_check/ceu_check_cxx_stdlib.hh"

#include <iostream>

int main()
{
    std::cout << ceu_interpret_cxx_std_version();
    std::cout << ceu_interpret_cxx_stdlib();
    std::cout << ceu_check_get_compiler_info();
    std::cout << ceu_check_get_ctypes_limit_info();
    std::cout << ceu_check_get_compile_time_os_info();
    std::cout << ceu_check_get_run_time_os_info();
}
