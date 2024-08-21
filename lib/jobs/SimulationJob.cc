#include "SimulationJob.hh"

namespace labw {
namespace art_modern {
    SimulationJob::~SimulationJob()
    {
        // delete fasta_fetch; // FIXME
    }

    SimulationJob::SimulationJob(BaseFastaFetch* fasta_fetch, const CoverageInfo& coverage_info, const int job_id)
        : fasta_fetch(fasta_fetch)
        , coverage_info(coverage_info)
        , job_id(job_id)
    {
    }
}
}