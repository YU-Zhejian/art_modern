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
