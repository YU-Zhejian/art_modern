#include "libam/utils/dump_utils.hh"

#include <boost/filesystem/operations.hpp>
#include <boost/stacktrace/safe_dump_to.hpp>
#include <boost/stacktrace/stacktrace.hpp>

#include <csignal>
#include <fstream>
#include <iostream>

namespace labw::art_modern {
constexpr char DUMP_FILENAME[] = "./backtrace.dump";
namespace {

    void my_signal_handler(const int signum)
    {
        std::signal(signum, SIG_DFL);
        boost::stacktrace::safe_dump_to(DUMP_FILENAME);
        std::raise(SIGABRT);
    }

} // namespace

void handle_dumps()
{
    std::signal(SIGSEGV, &my_signal_handler);
    std::signal(SIGABRT, &my_signal_handler);
    if (boost::filesystem::exists(DUMP_FILENAME)) {
        // there is a backtrace
        std::ifstream ifs(DUMP_FILENAME);

        const boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(ifs);
        std::cout << "Previous run crashed:\n" << st << std::endl;

        // cleaning up
        ifs.close();
        boost::filesystem::remove(DUMP_FILENAME);
    }
}
} // namespace labw::art_modern