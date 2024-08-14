#pragma once
#include <exception>

#define UNKNOWN_C_EXCEPTION "UNKNOWN"

namespace labw {
namespace art_modern {

    class CExceptionsProxy : public std::exception {
    public:
        enum class EXPECTATION { ZERO, NON_NEGATIVE, POSITIVE };

        explicit CExceptionsProxy(std::string c_lib_name, std::string details);
        const char* what() const noexcept override;

        static int requires_numeric(int c_value, std::string c_lib_name = UNKNOWN_C_EXCEPTION,
            const std::string& details = UNKNOWN_C_EXCEPTION, bool explain_using_strerror = false,
            EXPECTATION expectation = EXPECTATION::ZERO);
        static void* requires_not_null(void* c_value, std::string c_lib_name = UNKNOWN_C_EXCEPTION,
            const std::string& details = UNKNOWN_C_EXCEPTION, bool explain_using_strerror = false);

    private:
        std::string c_lib_name_ = UNKNOWN_C_EXCEPTION;
        std::string details_ = UNKNOWN_C_EXCEPTION;
    };

} // art_modern
} // labw
