#pragma once
#include "art_profile_builder/exe/APBConfig.hh"
#include "art_profile_builder/exe/IntermediateEmpDist.hh"

#include <memory>

namespace labw::art_modern {

void view_sam(const std::shared_ptr<IntermediateEmpDist>& ied1, const std::shared_ptr<IntermediateEmpDist>& ied2,
    const std::size_t thread_id, const APBConfig& config);
}
