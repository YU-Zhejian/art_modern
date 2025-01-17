#include "profile_utils.hh"

#include <ceu_check/ceu_check_os_macro.h>

#ifdef CEU_ON_POSIX
#include <sys/resource.h>
#endif

#include <boost/log/trivial.hpp>

namespace lab::art_modern::details {
void print_memory_usage(const char* file, const int line)
{
#ifdef CEU_ON_POSIX
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        BOOST_LOG_TRIVIAL(info) << file << ":" << line << ": " << "Memory usage: " << usage.ru_maxrss * 1024;
    }
#endif
}

} // namespace lab::art_modern::details
