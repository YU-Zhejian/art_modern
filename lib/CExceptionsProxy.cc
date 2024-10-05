#include <boost/log/trivial.hpp>

#include <boost/stacktrace.hpp>
#include <sstream>
#include <utility>

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
        return "Error occurred in C library";
    }

    CExceptionsProxy::CExceptionsProxy(std::string c_lib_name, std::string details)
        : c_lib_name_(std::move(c_lib_name))
        , details_(std::move(details))
    {
    }

} // art_modern
} // labw