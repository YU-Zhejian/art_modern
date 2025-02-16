#pragma once

#include "art_modern_config.h"

#ifdef WITH_BOOST_STACKTRACE
#include <boost/exception/all.hpp> // NOLINT
#include <boost/stacktrace/stacktrace.hpp>
#endif

namespace labw::art_modern {
template <class E> [[noreturn]] void throw_with_trace(const E& e)
{
#ifdef WITH_BOOST_STACKTRACE
    throw boost::enable_error_info(e) << boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>(
        boost::stacktrace::stacktrace());
#else
    throw e;
#endif
}
} // namespace labw::art_modern
