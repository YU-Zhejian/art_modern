/**
 * Copyright 2025-2026 YU Zhejian <yuzj25@seas.upenn.edu>
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
#include "libam_support/jobs/Scheduler.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <thread>

namespace labw::art_modern {

class JobPoolReporter : public Scheduler<std::chrono::seconds> {
public:
    JobPoolReporter(JobPool& jp, std::size_t reporting_interval_seconds);
    void callback() override;
    ~JobPoolReporter() override;
    DELETE_MOVE(JobPoolReporter)
    DELETE_COPY(JobPoolReporter)

private:
    JobPool& jp_;
    std::atomic<bool> should_stop_ { false };
    std::thread thread_;
};

} // namespace labw::art_modern
