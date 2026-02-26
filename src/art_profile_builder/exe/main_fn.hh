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
#include "art_modern_config.h" // NOLINT: For USED_HTSLIB_NAME

#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include <cstdlib>
#include <memory>
#include <vector>

namespace labw::art_modern {
void view_sam_mt(const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied1s,
    const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied2s, std::size_t n_threads, const APBConfig& config);

#ifdef WITH_NCBI_NGS
void view_sra_mt(const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied1s,
    const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied2s, std::size_t n_threads, const APBConfig& config);
#endif

} // namespace labw::art_modern
