#include "libam_support/jobs/JobPoolReporter.hh"

#include "libam_support/jobs/JobPool.hh"

#include <boost/log/trivial.hpp>

#include <cstdlib>
#include <thread>

namespace labw::art_modern {
JobPoolReporter::JobPoolReporter(labw::art_modern::JobPool &jp, const std::size_t reporting_interval_seconds)
        : jp_(jp)
        , reporting_interval_seconds_(reporting_interval_seconds)
{
}

        void JobPoolReporter::start() { thread_ = std::thread(&JobPoolReporter::job_, this); }
void JobPoolReporter::stop()
{
    should_stop_ = true;
    thread_.join();
}
void JobPoolReporter::job_() const

{
    std::this_thread::sleep_for(std::chrono::seconds(reporting_interval_seconds_));
    while (!should_stop_) {
        BOOST_LOG_TRIVIAL(info) << "JobPoolReporter: " << jp_.n_running_executors() << " JobExecutors running";
        std::this_thread::sleep_for(std::chrono::seconds(reporting_interval_seconds_));
    }
}
    } // namespace labw::art_modern