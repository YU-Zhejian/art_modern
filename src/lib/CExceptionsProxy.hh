/**
 *  @brief  C exceptions proxy.
 *  An exception proxy for C exceptions, with helper functions that asserts return value of C routines.
 */
#pragma once
#include "utils/exception_utils.hh"
#include <cerrno>
#include <cstring>
#include <exception>
#include <sstream>
#include <string>

namespace labw::art_modern {
static const std::string UNKNOWN_C_EXCEPTION = "UNKNOWN";

/*!
 * An exception proxy for C exceptions.
 */
class CExceptionsProxy : public std::exception {
public:
    /** Expected type of return value of C routines **/
    enum class EXPECTATION {
        /** Return value should be zero **/
        ZERO,
        /** Return value should be non-negative **/
        NON_NEGATIVE,
        /** Return value should be positive **/
        POSITIVE
    };

    CExceptionsProxy(std::string c_lib_name, std::string details);
    [[nodiscard]] const char* what() const noexcept override;
    void log() const;

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
    static t assert_numeric(const t c_value, const std::string& c_lib_name = UNKNOWN_C_EXCEPTION,
        const std::string& details = UNKNOWN_C_EXCEPTION, const bool explain_using_strerror = false,
        const EXPECTATION expectation = EXPECTATION::ZERO, const bool log = true)
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
            auto cep = CExceptionsProxy(c_lib_name, oss.str());
            if (log) {
                cep.log();
            }
            throw_with_trace(cep);
        }
        return c_value;
    }

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
    static t assert_not_null(const t c_value, const std::string& c_lib_name = UNKNOWN_C_EXCEPTION,
        const std::string& details = UNKNOWN_C_EXCEPTION, const bool explain_using_strerror = false,
        const bool log = true)
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

            auto cep = CExceptionsProxy(c_lib_name, oss.str());
            if (log) {
                cep.log();
            }
            throw_with_trace(cep);
        }
        return c_value;
    }

private:
    std::string c_lib_name_ = UNKNOWN_C_EXCEPTION;
    std::string details_ = UNKNOWN_C_EXCEPTION;
};

} // namespace labw::art_modern
