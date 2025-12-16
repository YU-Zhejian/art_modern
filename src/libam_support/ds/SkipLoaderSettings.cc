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

#include "libam_support/ds/SkipLoaderSettings.hh"

#include "libam_support/utils/mpi_utils.hh"

#include <cstdlib>

namespace labw::art_modern {
SkipLoaderSettings SkipLoaderSettings::from_mpi()
{
    if (have_mpi()) {
        return SkipLoaderSettings(mpi_size(), mpi_rank());
    }
    return SkipLoaderSettings(1, 0);
}

SkipLoaderSettings::SkipLoaderSettings(const std::size_t num_parallel_jobs, const std::size_t job_id)
    : num_parallel_jobs_(num_parallel_jobs)
    , job_id_(job_id)
{
}

std::size_t SkipLoaderSettings::skip_first() const { return job_id_; }

std::size_t SkipLoaderSettings::skip_others() const { return num_parallel_jobs_ - 1; }

} // namespace labw::art_modern
