#pragma once
#include <boost/exception/all.hpp>
#include <boost/stacktrace.hpp>

namespace labw {
namespace art_modern {

    using traced = boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>;

    template <class E> [[noreturn]] void throw_with_trace(const E& e)
    {
        throw boost::enable_error_info(e) << traced(boost::stacktrace::stacktrace());
    }
} // art_modern
} // labw
