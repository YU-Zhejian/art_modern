#pragma once
#include "ExceptionUtils.hh"
#include <exception>
#include <string>

namespace labw::art_modern {
const char UNKNOWN_C_EXCEPTION[] = "UNKNOWN";

class CExceptionsProxy : public std::exception {
public:
    enum class EXPECTATION { ZERO, NON_NEGATIVE, POSITIVE };

    CExceptionsProxy(std::string c_lib_name, std::string details);
    [[nodiscard]] const char* what() const noexcept override;
    void log() const;

    template <typename t>
    static t assert_numeric(t c_value, std::string c_lib_name = UNKNOWN_C_EXCEPTION,
        const std::string& details = UNKNOWN_C_EXCEPTION, bool explain_using_strerror = false,
        EXPECTATION expectation = EXPECTATION::ZERO, bool log = true)
    {
        if ((expectation == EXPECTATION::ZERO && c_value != 0) || (expectation == EXPECTATION::POSITIVE && c_value <= 0)
            || (expectation == EXPECTATION::NON_NEGATIVE && c_value < 0)) {
            std::ostringstream oss;
            if (details != UNKNOWN_C_EXCEPTION) {
                oss << details;
            }
            if (explain_using_strerror) {
                oss << strerror(errno);
            }
            oss << " returned " << c_value;
            auto cep = CExceptionsProxy(std::move(c_lib_name), oss.str());
            if (log) {
                cep.log();
            }
            throw_with_trace(cep);
        }
        return c_value;
    }

    template <typename t>
    static t assert_not_null(t c_value, std::string c_lib_name = UNKNOWN_C_EXCEPTION,
        const std::string& details = UNKNOWN_C_EXCEPTION, bool explain_using_strerror = false, bool log = true)
    {
        if (c_value == nullptr) {
            std::ostringstream oss;
            if (details != UNKNOWN_C_EXCEPTION) {
                oss << details;
            }
            if (explain_using_strerror) {
                oss << strerror(errno);
            }
            oss << " returned null";

            auto cep = CExceptionsProxy(std::move(c_lib_name), oss.str());
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
