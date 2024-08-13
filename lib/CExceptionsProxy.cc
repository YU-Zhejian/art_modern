#include <boost/log/trivial.hpp>

#include <boost/stacktrace.hpp>
#include <sstream>
#include <utility>

#include <cerrno>
#include <cstring>

#include "CExceptionsProxy.hh"

namespace labw {
namespace art_modern {
    const char* CExceptionsProxy::what() const noexcept
    {
        std::ostringstream oss;
        oss << "Error occurred in C library '" << c_lib_name_ << "' due to '" << details_ << "'";
        auto err_str = oss.str();
        BOOST_LOG_TRIVIAL(fatal) << err_str;
        BOOST_LOG_TRIVIAL(fatal) << "Stack trace:";
        BOOST_LOG_TRIVIAL(fatal) << boost::stacktrace::stacktrace();
        return err_str.c_str();
    }

    CExceptionsProxy::CExceptionsProxy(std::string c_lib_name, std::string details)
        : c_lib_name_(std::move(c_lib_name))
        , details_(std::move(details))
    {
    }
    int CExceptionsProxy::requires_numeric(int c_value,
        std::string c_lib_name,
        const std::string& details,
        bool explain_using_strerror, EXPECTATION expectation)
    {
        if ((expectation == EXPECTATION::ZERO && c_value != 0) || (expectation == EXPECTATION::POSITIVE && c_value <= 0) || (expectation == EXPECTATION::NON_NEGATIVE && c_value < 0)) {
            std::ostringstream oss;
            if (details == UNKNOWN_C_EXCEPTION) {
                oss << details;
            } else {
                if (explain_using_strerror) {
                    oss << strerror(errno);
                } else {
                    oss << "returned " << c_value;
                }
            }
            throw CExceptionsProxy(std::move(c_lib_name), oss.str());
        }
        return c_value;
    }

    void* CExceptionsProxy::requires_not_null(void* c_value, // NOSONAR cpp:S5008
        std::string c_lib_name,
        const std::string& details,
        bool explain_using_strerror)
    {
        if (c_value == nullptr) {
            std::ostringstream oss;
            if (details == UNKNOWN_C_EXCEPTION) {
                oss << details;
            } else {
                if (explain_using_strerror) {
                    oss << strerror(errno);
                } else {
                    oss << "returned null";
                }
            }
            throw CExceptionsProxy(std::move(c_lib_name), oss.str());
        }
        return c_value;
    }

} // art_modern
} // labw