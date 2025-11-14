#include "libam_support/jobs/JobPoolReporter.hh"

#include "libam_support/jobs/JobPool.hh"

#include <boost/log/trivial.hpp>

#include <chrono>
#include <cstdlib>
#include <thread>

namespace labw::art_modern {
JobPoolReporter::JobPoolReporter(labw::art_modern::JobPool& jp, const std::size_t reporting_interval_seconds)
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
    std::size_t current_sleep = 0;
    while (!should_stop_ && current_sleep < reporting_interval_seconds_) {
        current_sleep++;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    while (!should_stop_) {
        current_sleep = 0;
        BOOST_LOG_TRIVIAL(info) << "JobPoolReporter: " << jp_.n_running_executors() << " JobExecutors running";
        while (!should_stop_ && current_sleep < reporting_interval_seconds_) {
            current_sleep++;
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }
}
} // namespace labw::art_modern
