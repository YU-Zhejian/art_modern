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
#include "libam_support/jobs/JobPoolReporter.hh"

#include "libam_support/jobs/JobPool.hh"
#include "libam_support/jobs/Scheduler.hh"

#include <boost/log/trivial.hpp>

#include <chrono>
#include <cstdlib>

namespace labw::art_modern {
JobPoolReporter::JobPoolReporter(JobPool& jp, std::size_t reporting_interval_seconds)
    : Scheduler(std::chrono::seconds(reporting_interval_seconds), std::chrono::seconds(reporting_interval_seconds))
    , jp_(jp)
{
}

void JobPoolReporter::callback()

{
    BOOST_LOG_TRIVIAL(info) << "JobPoolReporter: " << jp_.n_running_executors() << " JobExecutors running";
}

JobPoolReporter::~JobPoolReporter() { stop(); }
} // namespace labw::art_modern
