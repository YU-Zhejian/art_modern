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

#include "libam_support/ds/CoverageInfo.hh"
#include "libam_support/ref/fetch/BaseFastaFetch.hh"
#include "libam_support/utils/class_macros_utils.hh"

#include <cstdlib>
#include <memory>

namespace labw::art_modern {

/**
 * A generic simulation job containing necessary data.
 *
 * The simulation job is move-constructible only.
 */
class SimulationJob {
public:
    /** Constructor that copies shared pointers.
     * @param fasta_fetch Fasta fetcher.
     * @param coverage_info Coverage information.
     * @param job_id Job ID.
     */
    SimulationJob(const std::shared_ptr<BaseFastaFetch>& fasta_fetch,
        const std::shared_ptr<CoverageInfo>& coverage_info, std::size_t job_id);
    /** Constructor that moves shared pointers.
     * @param fasta_fetch Fasta fetcher.
     * @param coverage_info Coverage information.
     * @param job_id Job ID.
     */
    SimulationJob(std::shared_ptr<BaseFastaFetch>&& fasta_fetch, std::shared_ptr<CoverageInfo>&& coverage_info,
        std::size_t job_id);

    SimulationJob(SimulationJob&& other) noexcept;
    DELETE_MOVE_ASSIGNMENT(SimulationJob)
    DELETE_COPY(SimulationJob)
    DEFAULT_DESTRUCTOR(SimulationJob)

    std::shared_ptr<BaseFastaFetch> fasta_fetch;
    std::shared_ptr<CoverageInfo> coverage_info;
    const std::size_t job_id;
};

} // namespace labw::art_modern
