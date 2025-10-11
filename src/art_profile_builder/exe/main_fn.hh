#pragma once
#include "art_profile_builder/lib/APBConfig.hh"
#include "art_profile_builder/lib/IntermediateEmpDist.hh"

#include <memory>

namespace labw::art_modern {

void view_sam(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
    const std::size_t thread_id, const APBConfig& config);
}
