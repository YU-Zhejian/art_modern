#include "art_modern_config.h"
#include <boost/log/attributes/clock.hpp>
#include <boost/log/attributes/current_thread_id.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>

#include <iostream>
#define DUMP_FILENAME "./backtrace.dump"

#ifdef WITH_MPI
#include "utils/mpi_log.hh"
#endif
#include "log_utils.hh"

namespace logging = boost::log;

namespace labw::art_modern {

void init_logger()
{
    auto core = logging::core::get();
    core->remove_all_sinks();
    core->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
    core->add_global_attribute("ThreadID", boost::log::attributes::current_thread_id());
#ifndef CEU_CM_IS_DEBUG
    core->set_filter(logging::trivial::severity >= logging::trivial::info);
#endif
#ifdef WITH_MPI
    core->add_global_attribute("MPIRank", MPIRankLoggerAttribute());
    core->add_global_attribute("MPIHostName", MPIHostNameLoggerAttribute());
    auto sink = boost::log::add_console_log(std::cerr,
        boost::log::keywords::format
        = "[%TimeStamp%] [T=%ThreadID%@MPI=%MPIRank%:%MPIHostName%] %Severity%: %Message%");
#else
    auto sink = boost::log::add_console_log(
        std::cerr, boost::log::keywords::format = "[%TimeStamp%] [Tread=%ThreadID%] %Severity%: %Message%");
#endif
    core->add_sink(sink);
}
void init_file_logger(const std::string& log_dir) { }
} // art_modern
// labw