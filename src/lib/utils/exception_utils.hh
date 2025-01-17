#pragma once
#include <boost/exception/all.hpp> // NOLINT
#include <boost/stacktrace/stacktrace.hpp>

namespace labw::art_modern {
template <class E> [[noreturn]] void throw_with_trace(const E& e)
{
    throw boost::enable_error_info(e) << boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>(
        boost::stacktrace::stacktrace());
}
} // namespace labw::art_modern
