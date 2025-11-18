/**
 * Copyright 2025 YU Zhejian <yuzj25@seas.upenn.edu>
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org/licenses/>.
 **/
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
