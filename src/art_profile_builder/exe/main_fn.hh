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
