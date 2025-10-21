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
