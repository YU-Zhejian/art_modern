#pragma once
#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include <cstdlib>
#include <memory>
#include <vector>

namespace labw::art_modern {

void view_sam(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
    std::size_t thread_id, const APBConfig& config);
void view_sam_mt(const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied1s,
    const std::vector<std::shared_ptr<IntermediateEmpDist>>& ied2s, std::size_t n_threads, const APBConfig& config);
} // namespace labw::art_modern
