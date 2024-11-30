#include <boost/log/trivial.hpp>

#include <boost/stacktrace.hpp>
#include <sstream>
#include <utility>

#include "CExceptionsProxy.hh"

namespace labw::art_modern {
const char* CExceptionsProxy::what() const noexcept { return "Error occurred in C library"; }

CExceptionsProxy::CExceptionsProxy(std::string c_lib_name, std::string details)
    : c_lib_name_(std::move(c_lib_name))
    , details_(std::move(details))
{
}
void CExceptionsProxy::log() const
{
    BOOST_LOG_TRIVIAL(fatal) << "Error occurred in C library '" << c_lib_name_ << "' due to '" << details_ << "'";
    BOOST_LOG_TRIVIAL(fatal) << "Stack trace:";
    BOOST_LOG_TRIVIAL(fatal) << boost::stacktrace::stacktrace();
}

} // art_modern
// labw