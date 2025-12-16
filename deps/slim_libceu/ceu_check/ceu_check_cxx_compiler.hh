/*!
 * @file ceu_check_cxx_compiler.hh
 * @brief Header to get compiler information at compile time
 *
 * @see https://sourceforge.net/p/predef/wiki/Architectures/
 */
#pragma once

#include <string>

/*!
 * @brief Get compiler information.
 *
 * @return Returned buffer, should be freed manually.
 */
std::string ceu_check_get_cxx_compiler_info();
