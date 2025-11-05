/**
 *  @brief  C exceptions proxy.
 *  An exception proxy for C exceptions, with helper functions that asserts return value of C routines.
 *
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#pragma once
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/mpi_utils.hh"

#include <boost/log/trivial.hpp>

#include <cerrno>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <string>

namespace labw::art_modern {
static const std::string UNKNOWN_C_EXCEPTION = "UNKNOWN";

/**
 * An exception proxy for C exceptions.
 * The current implementation is stupid. It generates horrific logs.
 */
class CExceptionsProxy {
public:
    /** Expected type of return value of C routines **/
    enum class EXPECTATION : std::uint8_t {
        /** Return value should be zero **/
        ZERO,
        /** Return value should be non-negative **/
        NON_NEGATIVE,
        /** Return value should be positive **/
        POSITIVE
    };

    /**
     * Default Constructor.
     */
    CExceptionsProxy() = delete;
    DELETE_COPY(CExceptionsProxy)
    DELETE_MOVE(CExceptionsProxy)
    DELETE_DESTRUCTOR(CExceptionsProxy)

    /**
     * Assert a C routine return value.
     *
     * @tparam t Some type that can be converted to int.
     * @param c_value Value returned by the C routine.
     * @param c_lib_name Name of the C library.
     * @param details Details about the error.
     * @param explain_using_strerror If true, the error message will include the error string returned by strerror.
     * @param expectation The expected value of the C routine.
     * @param log If true, the exception will be logged.
     * @return The `c_value`.
     */
    template <typename t>
    static t assert_numeric(t c_value, const std::string& c_lib_name = UNKNOWN_C_EXCEPTION,
        const std::string& details = UNKNOWN_C_EXCEPTION, bool explain_using_strerror = false,
        EXPECTATION expectation = EXPECTATION::ZERO, bool log = true);

    /**
     * Assert a C routine return value.
     *
     * @tparam t Some type that can should not be `null`.
     * @param c_value Value returned by the C routine.
     * @param c_lib_name Name of the C library.
     * @param details Details about the error.
     * @param explain_using_strerror If true, the error message will include the error string returned by strerror.
     * @param log If true, the exception will be logged.
     * @return The `c_value`.
     */
    template <typename t>
    static t assert_not_null(t c_value, const std::string& c_lib_name = UNKNOWN_C_EXCEPTION,
        const std::string& details = UNKNOWN_C_EXCEPTION, bool explain_using_strerror = false, bool log = true);
};

template <typename t>
t CExceptionsProxy::assert_numeric(const t c_value, const std::string& c_lib_name, const std::string& details,
    const bool explain_using_strerror, const EXPECTATION expectation, const bool log)
{
    if ((expectation == EXPECTATION::ZERO && c_value != 0) || (expectation == EXPECTATION::POSITIVE && c_value <= 0)
        || (expectation == EXPECTATION::NON_NEGATIVE && c_value < 0)) {
        std::ostringstream oss;
        if (details != UNKNOWN_C_EXCEPTION) {
            oss << details;
        }
        if (explain_using_strerror) {
            oss << std::strerror(errno);
        }
        oss << " returned " << c_value;
        if (log) {
            BOOST_LOG_TRIVIAL(fatal) << "Error occurred in C library '" << c_lib_name << "' due to '" << details << "'";
        }
        abort_mpi();
    }
    return c_value;
}

template <typename t>
t CExceptionsProxy::assert_not_null(const t c_value, const std::string& c_lib_name, const std::string& details,
    const bool explain_using_strerror, const bool log)
{
    if (c_value == nullptr) {
        std::ostringstream oss;
        if (details != UNKNOWN_C_EXCEPTION) {
            oss << details;
        }
        if (explain_using_strerror) {
            oss << std::strerror(errno);
        }
        oss << " returned null";

        if (log) {
            BOOST_LOG_TRIVIAL(fatal) << "Error occurred in C library '" << c_lib_name << "' due to '" << details << "'";
        }
        abort_mpi();
    }
    return c_value;
}

} // namespace labw::art_modern
