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

#include "libam_support/jobs/SimulationJob.hh"

#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"

#include <cstdlib>
#include <memory>
#include <utility>

namespace labw::art_modern {

SimulationJob::SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
    const std::shared_ptr<CoverageInfo>& coverage_info, const std::size_t job_id)
    : fasta_fetch(fasta_fetch)
    , coverage_info(coverage_info)
    , job_id(job_id)
{
}
SimulationJob::SimulationJob(std::shared_ptr<BaseFastaFetch>&& fasta_fetch,
    std::shared_ptr<CoverageInfo>&& coverage_info, const std::size_t job_id)
    : fasta_fetch(std::move(fasta_fetch))
    , coverage_info(std::move(coverage_info))
    , job_id(job_id)
{
}
SimulationJob::SimulationJob(SimulationJob&& other) noexcept
    : fasta_fetch(std::move(other.fasta_fetch))
    , coverage_info(std::move(other.coverage_info))
    , job_id(other.job_id)
{
}
} // namespace labw::art_modern
