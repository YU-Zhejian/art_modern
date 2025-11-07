/**
 * Copyright 2024-2025 YU Zhejian <yuzj25@seas.upenn.edu>
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

#include "libam_support/utils/class_macros_utils.hh"

#include <string>

namespace labw::art_modern {
/**
 * Abstract base class for job executors.

*/
class JobExecutor {
public:
    /**
     * Default constructor.
     */
    JobExecutor() = default;
    /** Default destructor */
    virtual ~JobExecutor() = default;
    DELETE_COPY(JobExecutor)
    DELETE_MOVE(JobExecutor)

    /**
     * Implement the job execution logic here.
     */
    virtual void operator()() = 0;
    /**
     * Implement to check whether the job is still running.
     */
    [[nodiscard]] virtual bool is_running() const = 0;
    /**
     * Implement to provide thread info for logging.
     */
    [[nodiscard]] virtual std::string thread_info() const = 0;
};

} // namespace labw::art_modern
