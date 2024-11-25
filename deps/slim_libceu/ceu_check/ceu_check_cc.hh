/*!
 * @file ceu_check_cc.h
 * @brief Header to get compiler information at compile time
 *
 * @see https://sourceforge.net/p/predef/wiki/Architectures/
 */
#pragma once

#include <ceu_check/ceu_check_cc_macro.h>
#include <string>

/*!
 * @brief Get compiler information.
 *
 * @return Returned buffer, should be freed manually.
 */
std::string ceu_check_get_compiler_info(void);
