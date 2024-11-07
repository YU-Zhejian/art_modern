#include "art_modern_config.h"
#include <boost/log/attributes/clock.hpp>
#include <boost/log/attributes/current_thread_id.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <boost/filesystem.hpp>

#include <iostream>

#ifdef WITH_MPI
#include "utils/mpi_log_attributes.hh"
#include <mpi.h>
#endif
#include "log_utils.hh"

namespace logging = boost::log;
namespace expr = logging::expressions;

namespace labw::art_modern {

void init_logger()
{
    auto core = logging::core::get();
    core->remove_all_sinks();
    core->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
    core->add_global_attribute("ThreadID", boost::log::attributes::current_thread_id());
#ifdef WITH_MPI
    core->add_global_attribute("MPIRank", MPIRankLoggerAttribute());
    core->add_global_attribute("MPIHostName", MPIHostNameLoggerAttribute());
#endif
    auto sink
        = boost::log::add_console_log(std::cerr, boost::log::keywords::format = "[%TimeStamp%] %Severity%: %Message%",
            logging::keywords::filter
            = expr::attr<int>("MPIRank") == 0 & logging::trivial::severity >= logging::trivial::info);

    core->add_sink(sink);
}
void init_file_logger(const std::string& log_dir)
{
    if (!boost::filesystem::exists(log_dir)) {
        boost::filesystem::create_directories(log_dir);
    }
#ifdef WITH_MPI
    int mpi_num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
    for (int i = 0; i < mpi_num_procs; i++) {
        std::stringstream file_name_ss;
        file_name_ss << log_dir << "/" << i << ".log";
        logging::add_file_log(logging::keywords::file_name = file_name_ss.str(),
            logging::keywords::format
            = "[%TimeStamp%] [T=%ThreadID%@MPI=%MPIRank%:%MPIHostName%] %Severity%: %Message%",
            logging::keywords::filter = expr::attr<int>("MPIRank") == i);
    }
    std::stringstream file_name_ss;
    file_name_ss << log_dir << "/"
                 << "nompi"
                 << ".log";
    logging::add_file_log(logging::keywords::file_name = file_name_ss.str(),
        logging::keywords::format = "[%TimeStamp%] [T=%ThreadID%@MPI=%MPIRank%:%MPIHostName%] %Severity%: %Message%",
        logging::keywords::filter = expr::attr<int>("MPIRank") == -1);
#else
    std::stringstream file_name_ss;
    file_name_ss << log_dir << "/"
                 << "nompi"
                 << ".log";
    logging::add_file_log(logging::keywords::file_name = file_name_ss.str(),
        logging::keywords::format = "[%TimeStamp%] [T=%ThreadID%] %Severity%: %Message%",
        logging::keywords::filter = expr::attr<int>("MPIRank") == -1);
#endif
}
}