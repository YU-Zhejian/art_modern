#include "art_modern_config.h"

#include "libam_support/CExceptionsProxy.hh"

#include <boost/log/trivial.hpp>
#ifdef WITH_BOOST_STACKTRACE
#include <boost/stacktrace/stacktrace.hpp>
#endif

#include <string>
#include <utility>

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
#ifdef WITH_BOOST_STACKTRACE
    BOOST_LOG_TRIVIAL(fatal) << "Stack trace:";
    BOOST_LOG_TRIVIAL(fatal) << boost::stacktrace::stacktrace();
#endif
}

} // namespace labw::art_modern
