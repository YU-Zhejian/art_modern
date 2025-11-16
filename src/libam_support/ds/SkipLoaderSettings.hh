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

#include <cstdlib>

namespace labw::art_modern {
class SkipLoaderSettings {
public:
    explicit SkipLoaderSettings(std::size_t num_parallel_jobs = 1, std::size_t job_id = 0);

    static SkipLoaderSettings from_mpi();

    /**
     * No. elements to skip at the beginning.
     * @return
     */
    [[nodiscard]] std::size_t skip_first() const;
    /**
     * No. elements to skip at the every period.
     * @return
     */
    [[nodiscard]] std::size_t skip_others() const;

private:
    /**
     * For MPI: Word size.
     */
    std::size_t num_parallel_jobs_;
    /**
     * For MPI: The rank.
     */
    std::size_t job_id_;
};

} // namespace labw::art_modern
