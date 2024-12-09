#pragma once
#include <ceu_check/ceu_check_c_cxx_std_macro.h>
#include <string>

/*!
 * @brief Get a nicely formatted compile-time C standard version number.
 */
std::string ceu_interpret_c_std_version();
/*!
 * @brief Get a nicely formatted compile-time CXX standard version number.
 */
std::string ceu_interpret_cxx_std_version();
