#pragma once

#include "libam_support/jobs/JobPool.hh"

#include <atomic>
#include <cstdlib>
#include <thread>

namespace labw::art_modern {

class JobPoolReporter {
public:
    explicit JobPoolReporter(JobPool& jp, std::size_t reporting_interval_seconds);
    void stop();
    void start();

private:
    void job_() const;

    JobPool& jp_;
    std::atomic<bool> should_stop_ { false };
    std::thread thread_;
    const std::size_t reporting_interval_seconds_;
};

} // namespace labw::art_modern
