#include "dump_utils.hh"

#include <boost/filesystem.hpp>
#include <boost/stacktrace.hpp>
#include <csignal>
#include <fstream>
#include <iostream>
namespace labw::art_modern {
const char DUMP_FILENAME[] = "./backtrace.dump";
void my_signal_handler(int signum)
{
    ::signal(signum, SIG_DFL);
    boost::stacktrace::safe_dump_to(DUMP_FILENAME);
    ::raise(SIGABRT);
}

void handle_dumps()
{
    ::signal(SIGSEGV, &my_signal_handler);
    ::signal(SIGABRT, &my_signal_handler);
    if (boost::filesystem::exists(DUMP_FILENAME)) {
        // there is a backtrace
        std::ifstream ifs(DUMP_FILENAME);

        boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(ifs);
        std::cout << "Previous run crashed:\n" << st << std::endl;

        // cleaning up
        ifs.close();
        boost::filesystem::remove(DUMP_FILENAME);
    }
}
}